import random

def calc_factors(num):
    factors = {}

    while num != 1:
        one_factor = calc_one_factor(num)
        factors[one_factor] = factors.get(one_factor, 0) + 1
        num = num // one_factor
    return factors

def calc_one_factor(num):
    sq_root = int(num**(1/2) + 1)
    for i in range(2, sq_root):
        if (num % i) == 0:
            return i
    return num

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


def modular_exponent(a, b, n):
    mask = 1
    result = 1
    while mask <= b:
        if mask & b:
            result = (result * a) % n
        a = (a * a) % n
        mask = mask << 1
    return result


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

    def check_solve(self, x_0, d):
        for i in range(d):
            x = (x_0 + self.n // d * i) % self.n
            # 遍历d个可能的解
            if modular_exponent(self.alpha, x, self.p) == self.beta:
                return x
        raise ValueError("算法失败")

    def solve(self):
        inv = None
        while not inv:
            init_a = random.randint(0, self.n)
            init_b = random.randint(0, self.n)
            init_x = modular_exponent(
                self.alpha, init_a, self.p)*modular_exponent(self.beta, init_b, self.p) % self.p

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
                inv = calc_inverse(self.n, b_2i - b_i)
            except ValueError as v_error:
                (d, c_0) = v_error.args
                if d < 1000:
                    # d是gcd(b_2i - b_i, n)，c_0是方程c * (b_2i - b_i) + y*n = d的一个特解
                    return self.check_solve(c_0 * (a_i - a_2i) // d, d)
                    # 乘以(a_i - a_2i) // d获得方程c * (b_2i - b_i) + y*n = (a_i - a_2i)的一个特解
        return (a_i - a_2i)*inv % self.n


class Pohlig_Hellman_algorithm:
    def __init__(self, factors: dict, alpha, beta, n, p):
        self.factors = factors
        self.alpha = alpha
        self.beta = beta
        self.n = n
        self.p = p

    def solve_one_factor(self, q, c):
        j = 0
        beta_j = self.beta
        numbers = []
        while j <= c - 1:
            sigma = modular_exponent(beta_j, self.n // (q ** (j + 1)), self.p)
            alpha = modular_exponent(self.alpha, (self.n // q), self.p)

            if q > 1000:
                solver = pollard_algorithm(alpha, sigma, q, self.p)
                a_j = solver.solve()
            else:
                for i in range(q):
                    if modular_exponent(alpha, i, self.p) == sigma:
                        a_j = i
                        break

            alpha_inv = calc_inverse(self.p, self.alpha)
            beta_j = (beta_j * modular_exponent(alpha_inv, a_j*(q**j), self.p)) % self.p
            j += 1
            numbers.append(a_j)

        walker = 1
        ret = 0
        for num in numbers:
            ret += walker*num
            walker = walker * q

        return ret

    def solve(self):
        M_i = []
        y_i = []
        a_i = []
        for (factor, power) in self.factors.items():
            the_M_i = self.n // (factor ** power)
            M_i.append(the_M_i)
            y_i.append(calc_inverse(factor ** power, the_M_i))
            a_i.append(self.solve_one_factor(factor, power))

        result = 0
        for (a, M, y) in zip(a_i, M_i, y_i):
            result += a*M*y
            result %= self.n

        return result
