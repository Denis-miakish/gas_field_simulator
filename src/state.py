from dataclasses import dataclass
from typing import Optional

@dataclass
class NodeState:
    
    p_in: float               # Давление на входе в элемент, атм
    p_out: float              # Давление на выходе из элемента, атм
    q_std: float              # Дебит газа при стандартных условиях, ст.м3/сут
    
    q_res: Optional[float] = None  # Дебит газа при пластовых условиях, м3/сут
    v: Optional[float] = None      # Средняя скорость газа, м/с
    rho: Optional[float] = None    # Средняя плотность газа, кг/м3
    mu: Optional[float] = None     # Динамическая вязкость газа, сП
    re: Optional[float] = None     # Число Рейнольдса (безразмерное)
