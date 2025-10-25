import math
from datetime import datetime

# Função original do Problema A
def f(x):
    return math.exp(-x) - 2 * math.sqrt(x)

# Forma contrativa escolhida (g1)
def g1(x):
    return 0.25 * math.exp(-2 * x)

# Método de Wegstein aplicado ao ponto fixo g(x)
def wegstein(g, f, x0, tol_abs=1e-6, tol_fun=1e-6, max_iter=100):
    print("============================================================")
    print("PROBLEMA A: f(x) = exp(-x) - 2*sqrt(x)")
    print("Método: Wegstein com Aceleração")
    print(f"Executado em: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("============================================================")
    print(f"x0 = {x0:.6f}")
    print(f"Tolerâncias: |Δx| ≤ {tol_abs}, |f(x)| ≤ {tol_fun}")
    print("------------------------------------------------------------")

    # Semente inicial coerente com ponto fixo
    x_prev = x0
    x_curr = g(x_prev)
    iter_count = 1

    # Iterações de Wegstein
    while iter_count < max_iter:
        g_prev = g(x_prev)
        g_curr = g(x_curr)

        # Evita divisão por zero
        if abs(x_curr - x_prev) < 1e-14:
            print("Divisão por zero detectada no cálculo de s_k. Encerrando.")
            break

        # Cálculo de s_k e w_k (fatores de aceleração)
        s_k = (g_curr - g_prev) / (x_curr - x_prev)
        w_k = 1.0 / (1.0 - s_k)

        # Limite de estabilidade (clamp do fator de aceleração)
        w_k = max(min(w_k, 1.2), 0.2)

        # Próxima iteração com aceleração
        x_next = (1 - w_k) * x_curr + w_k * g_curr

        # Projeção de domínio (garante x ≥ 0)
        if x_next < 0:
            x_next = abs(x_next) / 2

        # Critérios de parada
        dx = abs(x_next - x_curr)
        fx = abs(f(x_next))

        print(f"Iter {iter_count:02d} | x = {x_next:.8f} | |Δx|={dx:.2e} | |f(x)|={fx:.2e} | w={w_k:.3f}")

        if dx <= tol_abs and fx <= tol_fun:
            print("\nConvergência alcançada com sucesso!")
            return x_next, iter_count, fx, True

        # Atualiza valores
        x_prev, x_curr = x_curr, x_next
        iter_count += 1

    print("\nNão convergiu dentro do número máximo de iterações.")
    return x_curr, iter_count, abs(f(x_curr)), False


# Execução principal
if __name__ == "__main__":
    x0 = 0.2  # estimativa inicial coerente
    raiz, iters, fx, convergiu = wegstein(g1, f, x0)

    print("\n============================================================")
    print("RESULTADOS - PROBLEMA A (Método de Wegstein)")
    print("============================================================")
    print(f"Raiz estimada     : {raiz:.6f}")
    print(f"|f(raiz)|         : {fx:.3e}")
    print(f"Iterações         : {iters}")
    print(f"Status            : {'✓ Convergiu' if convergiu else '✗ Não convergiu'}")
    print("============================================================")
