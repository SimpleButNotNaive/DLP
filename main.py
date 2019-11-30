# import time
import args
import func

if __name__ == "__main__":
    ya_solver = func.Pohlig_Hellman_algorithm(func.calc_factors(
        args.ord_g), args.g, args.ya, args.ord_g, args.p)
    yb_solver = func.Pohlig_Hellman_algorithm(func.calc_factors(
        args.ord_g), args.g, args.yb, args.ord_g, args.p)

    ya_solver_pollard = func.pollard_algorithm(args.g, args.ya, args.ord_g, args.p)
    yb_solver_pollard = func.pollard_algorithm(args.g, args.yb, args.ord_g, args.p)


    xa = ya_solver.solve()
    xb = yb_solver.solve()

    print("ya的离散对数xa为:", xa)
    print("yb的离散对数xb为:", xb)

    print("ya^xb = \n", func.modular_exponent(args.ya, xb, args.p))
    print("yb^xa = \n", func.modular_exponent(args.yb, xa, args.p))
    print("两者相等，它们即DH密钥交换协议中的共同内容")

    # start = time.time()
    # for i in range(20):
    #     ya_solver.solve()
    # end = time.time()
    # print("使用Pohlig Hellman算法的运行时间：（计算十次）")
    # print(end - start)

    # start = time.time()
    # for i in range(20):
    #     ya_solver_pollard.solve()
    # end = time.time()
    # print("使用Pollard rho算法的运行时间：（计算十次）")
    # print(end - start)
