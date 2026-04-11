class LinearInterpolator:
    def __init__(self, xs: list[float], ys: list[float]):
        if len(xs) != len(ys):
            raise ValueError("Списки xs и ys должны быть одинаковой длины")
        self.xs = xs
        self.ys = ys

    def interpolate(self, xp: float) -> float:
        
        if xp < self.xs[0] or xp > self.xs[-1]:
            raise ValueError

        for i in range(len(self.xs) - 1):
            if self.xs[i] <= xp <= self.xs[i+1]:
                x_i, x_next = self.xs[i], self.xs[i+1]
                y_i, y_next = self.ys[i], self.ys[i+1]
                return y_i + (y_next - y_i) / (x_next - x_i) * (xp - x_i)

        return self.ys[-1]
