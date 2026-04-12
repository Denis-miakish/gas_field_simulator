from src.state import NodeState

class DCS:
    """
    Класс дожимной компрессорной станции (ДКС).
    Обеспечивает повышение давления газа до уровня магистрали.
    """
    def __init__(self, p_line: float, cr: float, q_ext: float):
        """
        :param p_line: Давление в магистрали (константа), атм [3, 4]
        :param cr: Степень сжатия (Compression Ratio) [3, 4]
        :param q_ext: Дебит стороннего газа, поступающего на манифолд, ст.м3/сут [3, 4]
        """
        self.p_line = p_line
        self.cr = cr
        self.q_ext = q_ext

    def get_p_in(self) -> float:
        """
        Расчет давления на входе в ДКС.
        P_in,DCS = P_line / CR [5, 6]
        """
        # Если CR = 1, станция отключена и давление на входе равно давлению в линии [7, 8]
        if self.cr <= 1:
            return self.p_line
        
        return self.p_line / self.cr

    def get_state(self, q_cluster_std: float) -> NodeState:
        """
        Формирует объект состояния для ДКС.
        :param q_cluster_std: Суммарный дебит трех скважин куста, ст.м3/сут
        """
        p_in = self.get_p_in()
        
        # Общий дебит через ДКС включает газ куста и сторонний газ [9, 10]
        q_total = q_cluster_std + self.q_ext
        
        # Согласно ТЗ, для ДКС поля q_res, v, rho могут быть None [11, 12]
        return NodeState(
            p_in=p_in,
            p_out=self.p_line,
            q_std=q_total,
            q_res=None,
            v=None,
            rho=None,
            re=None
        )
