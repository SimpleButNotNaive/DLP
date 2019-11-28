import random
import args


class pollard_algorithm:
    def __init__(self, alpha, beta, n, p):
        self.alpha = alpha
        self.beta = beta
        self.n = n
        self.p = p
        self.func_list = {1: self.func_1, 0: self.func_2, 2: self.func_3}

    def func_1(self, x, a, b):
        return (self.beta*x % self.p, a, (b + 1) % self.n)

    def func_2(self, x, a, b):
        return (x**2 % self.p, 2*a % self.n, 2*b % self.n)

    def func_3(self, x, a, b):
        return (self.alpha*x % self.p, (a + 1) % self.n, b)

    def func(self, x, a, b):
        return self.func_list[x % 3](x, a, b)

    @staticmethod
    def calc_inverse(n, ele):
        a = n
        b = ele % n
        t_0 = 0
        t = 1
        q = a // b
        r = a % b
        while r > 0:
            temp = (t_0 - q*t) % n
            t_0 = t
            t = temp
            a = b
            b = r
            q = a // b
            r = a % b
        if b != 1:
            raise ValueError(b, t)
            # 如果算法失败则有b = t*ele + kn
            # b是最大公因数，t是方程的一个特解
        return t

    @staticmethod
    def modular_exponent(a, b, n):
        mask = 1
        result = 1
        while mask <= b:
            if mask & b:
                result = (result * a) % n
            a = (a * a) % n
            mask = mask << 1
        return result

    def check_solve(self, x_0, d):
        for i in range(d):
            x = (x_0 + self.n // d * i) % self.n
            # 遍历d个可能的解
            if pollard_algorithm.modular_exponent(self.alpha, x, self.p) == self.beta:
                return x
        raise ValueError("算法失败")

    def solve(self):
        inv = None
        while not inv:
            init_a = random.randint(0, self.n)
            init_b = random.randint(0, self.n)
            init_x = pollard_algorithm.modular_exponent(
                self.alpha, init_a, self.p)*pollard_algorithm.modular_exponent(self.beta, init_b, self.p) % self.p

            tuple_1 = self.func(init_x, init_a, init_b)
            tuple_2 = self.func(*tuple_1)

            while tuple_1[0] != tuple_2[0]:
                tuple_1 = self.func(*tuple_1)
                tuple_2 = self.func(*tuple_2)
                tuple_2 = self.func(*tuple_2)

            a_i = tuple_1[1]
            a_2i = tuple_2[1]

            b_i = tuple_1[2]
            b_2i = tuple_2[2]

            try:
                inv = pollard_algorithm.calc_inverse(self.n, b_2i - b_i)
            except ValueError as v_error:
                (d, c_0) = v_error.args
                if d < 1000:
                    # d是gcd(b_2i - b_i, n)，c_0是方程c * (b_2i - b_i) + y*n = d的一个特解
                    return self.check_solve(c_0 * (a_i - a_2i) // d, d)
                    # 乘以(a_i - a_2i) // d获得方程c * (b_2i - b_i) + y*n = (a_i - a_2i)的一个特解
        return (a_i - a_2i)*inv % self.n


if __name__ == "__main__":
    p = pollard_algorithm(args.g, args.ya, args.ord_g, args.p)
    solution = p.solve()
    print(pollard_algorithm.modular_exponent(args.g, solution, args.p))
    print(args.ya)

    p = pollard_algorithm(args.g, args.yb, args.ord_g, args.p)
    solution2 = p.solve()
    print(pollard_algorithm.modular_exponent(args.g, solution2, args.p))
    print(args.yb)

    print(pollard_algorithm.modular_exponent(args.ya, solution2, args.p))
    print(pollard_algorithm.modular_exponent(args.yb, solution, args.p))
