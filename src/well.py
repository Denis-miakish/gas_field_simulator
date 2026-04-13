import math
from src.fluid import Fluid

class Well:
    """
    Класс скважины, описывающий приток газа из пласта (закон Дарси).
    """
    def __init__(self, k_md: float, h_m: float, re_m: float, rw_m: float, fluid: Fluid):
        """
        :param k_md: Проницаемость, мД
        :param h_m: Эффективная мощность пласта, м
        :param re_m: Радиус контура питания, м
        :param rw_m: Радиус скважины, м
        :param fluid: Объект класса Fluid для получения свойств газа
        """
        self.k = k_md
        self.h = h_m
        self.re = re_m
        self.rw = rw_m
        self.fluid = fluid
        
        # Коэффициент перевода единиц из ТЗ (раздел 2.3)
        self.beta = 0.00852702 

    def get_c(self, p_res: float) -> float:
        """
        Расчет коэффициента продуктивности C.
        C = (beta * k * h) / (mu * ln(re/rw))
        """
        mu = self.fluid.get_mu(p_res) # вязкость при пластовом давлении
        ln_term = math.log(self.re / self.rw)
        
        return (self.beta * self.k * self.h) / (mu * ln_term)

    def q_std(self, p_res: float, p_bhp: float) -> float:
        """
        Расчет дебита газа по закону Дарси (ст.м3/сут).
        q = C * (P_res - P_bhp)
        """
        if p_bhp >= p_res:
            return 0.0
            
        c = self.get_c(p_res)
        return c * (p_res - p_bhp)
