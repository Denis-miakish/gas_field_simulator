import numpy as np
from scipy.optimize import fsolve
from src.state import NodeState
from src.reservoir import Reservoir
from src.well import Well
from src.pipe import Pipe
from src.compressor import DCS


class FieldSimulator:
    """
    Класс-симулятор для нахождения гидравлического равновесия системы
    и моделирования разработки во времени.
    """

    def __init__(
        self,
        reservoir: Reservoir,
        wells: list[Well],
        well_pipes: list[Pipe],
        shlyf: Pipe,
        dcs: DCS,
    ):
        self.reservoir = reservoir
        self.wells = wells
        self.well_pipes = well_pipes
        self.shlyf = shlyf
        self.dcs = dcs

    def solve(self) -> dict[str, NodeState]:
        """
        Нахождение рабочей точки: q1, q2, q3, P_man.
        Возвращает словарь NodeState для well_1..3, shlyf, dcs.
        """
        p_res = self.reservoir.props.p  # текущее пластовое давление

        def equations(vars_):
            # x = [q1, q2, q3, P_man]
            q1, q2, q3, p_man = vars_
            qs = [q1, q2, q3]
            res = []

            # 1–3. Уравнения притока по Дарси для каждой скважины
            for i in range(3):
                # состояние НКТ при заданном q и P_man
                state_well = self.well_pipes[i].dp(qs[i], p_man)
                # перепад в трубе
                dp_tube = state_well.p_in - state_well.p_out
                p_bhp = p_man + dp_tube

                q_calc = self.wells[i].q_std(p_res, p_bhp)
                res.append(qs[i] - q_calc)

            # 4. Баланс давлений на манифольде
            q_total_shlyf = sum(qs) + self.dcs.q_ext
            p_in_dcs = self.dcs.get_p_in()

            state_shlyf = self.shlyf.dp(q_total_shlyf, p_in_dcs)
            dp_shlyf = state_shlyf.p_in - state_shlyf.p_out
            p_man_calc = p_in_dcs + dp_shlyf

            res.append(p_man - p_man_calc)
            return res

        # начальное приближение
        initial_guess = [
            500.0,
            500.0,
            500.0,
            self.dcs.get_p_in() + 5.0,
        ]

        sol = fsolve(equations, initial_guess)  # [q1, q2, q3, P_man]

        # дебиты неотрицательные
        q_final = [max(0.0, val) for val in sol[:3]]
        p_man_final = sol[3]  # ИСПРАВЛЕНО: правильный индекс

        p_in_dcs = self.dcs.get_p_in()

        states: dict[str, NodeState] = {}

        # состояния скважин: считаем через их трубы при найденных q и P_man
        for i in range(3):
            states[f"well_{i+1}"] = self.well_pipes[i].dp(q_final[i], p_man_final)

        # состояние шлейфа
        states["shlyf"] = self.shlyf.dp(sum(q_final) + self.dcs.q_ext, p_in_dcs)

        # состояние ДКС (как ты реализовал в DCS.get_state)
        states["dcs"] = self.dcs.get_state(sum(q_final))

        return states

    def run(self, n_days: int, dt: float = 1.0) -> list[dict]:
        """
        Динамическая симуляция на n_days, явная схема по времени.
        Возвращает список словарей с историей.
        """
        history: list[dict] = []

        for day in range(n_days):
            # 1. Решаем систему на текущем шаге
            step_state = self.solve()
            q_total_cluster = sum(step_state[f"well_{i+1}"].q_std for i in range(3))

            # 2. Сохраняем результаты
            gp_prev = history[-1]["gp"] if history else 0.0
            history.append(
                {
                    "day": day,
                    "p_res": self.reservoir.props.p,
                    "q_total": q_total_cluster,
                    "p_man": step_state["shlyf"].p_in,
                    "gp": gp_prev + q_total_cluster * dt,
                }
            )

            # 3. Обновляем пластовое давление (в баланс входит только добыча куста)
            self.reservoir.props.p = self.reservoir.p2(q_total_cluster, dt)

            if day % 10 == 0:
                print(f"Шаг {day}/{n_days}: P_res = {self.reservoir.props.p:.2f} атм")

        return history
