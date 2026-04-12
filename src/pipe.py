import math
from src.fluid import Fluid
from src.state import NodeState

class Pipe:
    """
    Класс трубопровода для расчета гидравлических потерь в НКТ и шлейфе.
    """
    def __init__(self, l_m: float, d_m: float, vertical_depth_m: float, 
                 roughness_m: float, fluid: Fluid):
        self.L = l_m                   # Длина, м
        self.D = d_m                   # Внутренний диаметр, м
        self.H = vertical_depth_m      # Вертикальная глубина, м (для шлейфа = 0)
        self.delta = roughness_m       # Шероховатость, м
        self.fluid = fluid
        self.g = 9.81                  # м/с2

    def _get_lambda(self, re: float) -> float:
        """Расчет коэффициента трения лямбда (Colebrook-White)."""
        if re < 1e-10: # Защита от нулевого дебита
            return 0.0
        
        # Ламинарный режим
        if re < 2300:
            return 64 / re
        
        # Турбулентный режим (итерации по Колбруку-Уайту)
        l_curr = 0.02 # Начальное приближение
        for _ in range(100):
            # Формула итерации из раздела 2.4 ТЗ
            term = self.delta / (3.7 * self.D) + 2.51 / (re * math.sqrt(l_curr))
            l_next = (-2 * math.log10(term))**(-2)
            
            if abs(l_next - l_curr) < 1e-6:
                return l_next
            l_curr = l_next
        return l_curr

    def dp(self, q_std: float, p_atm: float) -> NodeState:
        """
        Расчет перепада давления и состояния потока.
        :param q_std: Дебит при ст.у., ст.м3/сут
        :param p_atm: Давление в узле, атм
        :return: Объект NodeState
        """
        # 1. Получаем свойства газа при текущем давлении
        rho = self.fluid.get_rho(p_atm)
        mu_cp = self.fluid.get_mu(p_atm)
        mu_pas = mu_cp / 1000 # Перевод сП в Па*с для Re
        bg = self.fluid.get_bg(p_atm)
        
        # 2. Расчет скорости газа (раздел 2.4 ТЗ)
        # v = (4 * q_std * Bg) / (pi * D^2 * 86400)
        v = (4 * q_std * bg) / (math.pi * self.D**2 * 86400)
        
        # 3. Число Рейнольдса
        re = (rho * v * self.D) / mu_pas
        
        # 4. Коэффициент трения
        lam = self._get_lambda(re)
        
        # 5. Перепад давления по Дарси-Вейсбаху (в атм)
        # delta_p = (1 / 101325) * (friction + hydrostatic)
        friction_term = lam * (self.L / self.D) * (rho * v**2 / 2.0)
        hydrostatic_term = rho * self.g * self.H
        
        delta_p_pa = friction_term + hydrostatic_term
        delta_p_atm = delta_p_pa / 101325
        
        # 6. Формирование состояния (давление на выходе p_out условно)
        return NodeState(
            p_in=p_atm,
            p_out=p_atm - delta_p_atm,
            q_std=q_std,
            q_res=q_std * bg,
            v=v,
            rho=rho,
            mu=mu_cp,
            re=re
        )
