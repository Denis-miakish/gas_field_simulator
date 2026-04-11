import math
from src.interpolator import LinearInterpolator

class Fluid:
    def __init__(self, t_res: float, p_std_kgm3: float, xa_pct: float, xy_pct: float, 
                 visc_data_x: list[float], visc_data_y: list[float]):
       
        self.T = t_res
        self.rho_std = p_std_kgm3
        # Перевод процентов в доли 
        self.xa = xa_pct / 100
        self.xy = xy_pct / 100
        self.xe = 1.0 - self.xa - self.xy # Доля эквивалентного углеводорода 
        
        # Константы 
        self.R = 8.314
        self.P_std = 101325  # Па
        self.T_std = 293.15  # К
        
        # Расчет молярной массы всего газа 
        
        self.M = (self.rho_std * self.R * self.T_std) / self.P_std
        
        # Параметры для GERG-91 [17, 18]
        # Me - молярная масса эквивалентного углеводорода
        self.Me = (24.05525 * self.rho_std - 28.0135 * self.xa - 44.01 * self.xy) / self.xe
        self.H = 128.64 + 47.479 * self.Me # Теплота сгорания 
        
        # Инициализация интерполятора для вязкости 
        self.visc_interp = LinearInterpolator(visc_data_x, visc_data_y)

    def get_z(self, p_atm: float) -> float:
        """Расчет Z-фактора по GERG-91 (ГОСТ 30319.2-96) [15, 19]."""
        p_mpa = p_atm * 0.101325 # Перевод атм в МПа 
        T = self.T
        H = self.H
        
        
        # B-коэффициенты
        b1 = -0.425468 + 2.865e-3*T - 4.62073e-6*T**2 + \
             (8.77118e-4 - 5.56281e-6*T + 8.81514e-9*T**2)*H + \
             (-8.24747e-7 + 4.31436e-9*T - 6.08319e-12*T**2)*H**2
        b2 = -0.1446 + 7.4091e-4*T - 9.1195e-7*T**2
        b3 = -0.86834 + 4.0376e-3*T - 5.1657e-6*T**2
        b23 = -0.339693 + 1.61176e-3*T - 2.04429e-6*T**2
        b_star = 0.72 + 1.875e-5 * (320.0 - T)**2
        
        # С-коэффициенты
        c1 = -0.302488 + 1.95861e-3*T - 3.16302e-6*T**2 + \
             (6.46422e-4 - 4.22876e-6*T + 6.88157e-9*T**2)*H + \
             (-3.32805e-7 + 2.2316e-9*T - 3.67713e-12*T**2)*H**2
        c2 = 7.8498e-3 - 3.9895e-5*T + 6.1187e-8*T**2
        c3 = 2.0513e-3 + 3.4888e-5*T - 8.3703e-8*T**2
        c223 = 5.52066e-3 - 1.68609e-5*T + 1.57169e-8*T**2
        c233 = 3.58783e-3 + 8.06674e-6*T - 3.25798e-8*T**2
        c_star = 0.92 + 0.0013 * (T - 270.0)

        # Смесевые коэффициенты Bm и Cm 
        bm = (self.xe**2 * b1) + (self.xe * self.xa * b_star * (b1 + b2)) - \
             (1.73 * self.xe * self.xy * math.sqrt(abs(b1 * b3))) + \
             (self.xa**2 * b2) + (2.0 * self.xa * self.xy * b23) + (self.xy**2 * b3)
        
        cm = (self.xe**3 * c1) + (self.xa**3 * c2) + (self.xy**3 * c3) + \
             (3 * self.xa * self.xy**2 * c233) + (3 * self.xa**2 * self.xy * c223)
        
        # Решение кубического уравнения 
        b_val = 1000 * p_mpa / (2.7715 * T)
        B0 = b_val * bm
        C0 = b_val**2 * cm
        
        A1 = 1 + B0
        A0 = 1 + 1.5 * (B0 + C0)
        
        discr = A0**2 - A1**3
        if discr < 0: return 1 # Защита от некорректных давлений
        
        A2 = (A0 + math.sqrt(discr))**(1/3)
        return (1 + A2 + A1/A2) / 3

    def get_rho(self, p_atm: float) -> float:
        """Плотность газа: rho = (P * M) / (Z * R * T) [1, 15, 19]."""
        z = self.get_z(p_atm)
        p_pa = p_atm * 101325
        return (p_pa * self.M) / (z * self.R * self.T)

    def get_mu(self, p_atm: float) -> float:
        """Динамическая вязкость через интерполятор [2, 19, 22]."""
        return self.visc_interp.interpolate(p_atm)

    def get_bg(self, p_atm: float) -> float:
        """Объемный коэффициент: Bg = (P_std * Z * T) / (P * T_std) [19, 22, 23]."""
        z = self.get_z(p_atm)
        p_pa = p_atm * 101325
        return (self.P_std * z * self.T) / (p_pa * self.T_std)
