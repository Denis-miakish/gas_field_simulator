from dataclasses import dataclass
from src.fluid import Fluid

@dataclass
class ResProps:
    """
    Класс для хранения текущих свойств пласта.
    """
    p: float        # Текущее пластовое давление, атм
    v_res: float    # Объем пласта, м3
    t_res: float    # Пластовая температура, К

class Reservoir:
    """
    Класс пласта («бак»), рассчитывающий материальный баланс.
    """
    def __init__(self, props: ResProps, fluid: Fluid):
        self.props = props
        self.fluid = fluid

    def p2(self, q_total_std: float, dt: float = 1) -> float:
        """
        Расчет нового пластового давления по формуле материального баланса.
        
        :param q_total_std: Суммарный дебит всех скважин куста, ст.м3/сут
        :param dt: Шаг по времени, сут (по умолчанию 1 день)
        :return: Новое пластовое давление P_res^(t+1), атм
        """
        p_current = self.props.p
        v_res = self.props.v_res
        
        # Получаем свойства газа при текущем давлении
        z = self.fluid.get_z(p_current)
        rho = self.fluid.get_rho(p_current)
        rho_std = self.fluid.rho_std
        
        # Формула материального баланса из раздела 2.3 ТЗ:
        pressure_drop = (z * rho_std / rho) * (q_total_std / v_res) * dt
        
        p_next = p_current - pressure_drop
        
        # Важно: согласно ТЗ (п. 3.5), метод p2 только возвращает значение,
        # но не меняет состояние self.props.p напрямую.
        return p_next
