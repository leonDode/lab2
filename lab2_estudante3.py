#!/usr/bin/env python3
"""
Lab 2 - Métodos de Raízes Não Lineares
Estudante 3: Implementação de Substituições Sucessivas e Wegstein

Este módulo implementa métodos para encontrar raízes de equações não lineares:
- Substituições Sucessivas (Ponto Fixo)
- Wegstein (aceleração)

Problemas implementados:
- Problema A: f(x) = exp(-x) - 2*sqrt(x) no intervalo [0,1]
- Problema B: Equação de escoamento turbulento em tubo liso
"""

import math
import time
import csv
from dataclasses import dataclass
from typing import List, Tuple, Dict, Callable, Optional
from pathlib import Path


@dataclass
class RootResult:
    """Resultado de um método de busca de raiz."""
    root: float
    fval: float
    iters: int
    converged: bool
    history: List[Tuple[int, float, float]]  # (iter, x, f(x))


def is_finite(x: float) -> bool:
    """Verifica se x é finito (não NaN nem infinito)."""
    return math.isfinite(x)


def clamp_domain(x: float, min_val: float = 1e-10) -> float:
    """Aplica clamp para garantir domínio válido (ex: x > 0 para log)."""
    return max(x, min_val)


def _fmt_num(x: float) -> str:
    if x is None or not math.isfinite(x):
        return f"{x}"
    return f"{x:12.6f}"

def _fmt_res(r: float) -> str:
    if r is None or not math.isfinite(r):
        return f"{r}"
    return f"{abs(r):12.2e}"

def format_table_row(method: str, root: float, fval: float, iters: int, converged: bool) -> str:
    status = "✓" if converged else "✗"
    root_s = _fmt_num(root)
    fval_s = _fmt_res(fval)
    return f"{method:15} | {root_s} | {fval_s} | {iters:4d} | {status}"


def save_history_csv(filename: str, history: List[Tuple[int, float, float]], problem: str, method: str):
    """Salva histórico de iterações em CSV."""
    Path("./out").mkdir(exist_ok=True)
    filepath = f"./out/{filename}"
    
    with open(filepath, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Iteracao', 'x', 'f(x)', 'Problema', 'Metodo'])
        for iter_num, x, fx in history:
            writer.writerow([iter_num, x, fx, problem, method])


def fixed_point(
    g: Callable[[float], float],
    x0: float,
    eps_abs_x: float = 1e-6,
    delta_f: float = 1e-6,
    max_iter: int = 500,
    f_original: Optional[Callable[[float], float]] = None
) -> RootResult:
    """
    Iteração de ponto fixo: x_{k+1} = g(x_k).
    Convergência somente se: |x_{k+1}-x_k| <= eps_abs_x E |f(x_{k+1})| <= delta_f.
    """
    if f_original is None:
        f_original = g

    history: List[Tuple[int, float, float]] = []
    x = x0

    try:
        fx = f_original(x)
    except Exception:
        return RootResult(root=x, fval=float("nan"), iters=0, converged=False, history=history)
    history.append((0, x, fx))

    for k in range(1, max_iter + 1):
        try:
            x_new = g(x)
            if not is_finite(x_new):
                return RootResult(root=x, fval=fx, iters=k-1, converged=False, history=history)

            fx_new = f_original(x_new)  # sempre residual da FUNÇÃO ORIGINAL
            history.append((k, x_new, fx_new))

            step_ok = abs(x_new - x) <= eps_abs_x
            res_ok  = abs(fx_new)     <= delta_f

            if step_ok and res_ok:
                return RootResult(root=x_new, fval=fx_new, iters=k, converged=True, history=history)

            x = x_new
            fx = fx_new

        except (ValueError, ZeroDivisionError, OverflowError):
            # Fora de domínio ou numérico ruim => falhou
            return RootResult(root=x, fval=fx, iters=k-1, converged=False, history=history)

    return RootResult(root=x, fval=fx, iters=max_iter, converged=False, history=history)


def wegstein(
    g: Callable[[float], float],
    x0: float,
    x1: float,
    eps_abs_x: float = 1e-6,
    delta_f: float = 1e-6,
    max_iter: int = 500,
    w_bounds: Tuple[float, float] = (-10.0, 10.0),
    f_original: Optional[Callable[[float], float]] = None
) -> RootResult:
    """
    Wegstein sobre a iteração de ponto fixo g(x). 
    Usa w = -dk / (dk - dk_1). 
    Fallback seguro: se |dk - dk_1| pequeno => w = 1 (iteração simples).
    Convergência somente se passo E residual passarem.
    """
    if f_original is None:
        f_original = g

    history: List[Tuple[int, float, float]] = []

    # Garantir semente: x1 = g(x0) é o padrão correto
    try:
        if not is_finite(x1):
            x1 = g(x0)
    except Exception:
        return RootResult(root=x0, fval=float("nan"), iters=0, converged=False, history=history)

    # Avalia pontos
    try:
        fx0 = f_original(x0)
        history.append((0, x0, fx0))
        fx1 = f_original(x1)
        history.append((1, x1, fx1))
    except (ValueError, ZeroDivisionError, OverflowError):
        return RootResult(root=x0, fval=float("nan"), iters=0, converged=False, history=history)

    x_prev, x_curr = x0, x1
    try:
        g_prev = g(x_prev)
        g_curr = g(x_curr)
    except Exception:
        return RootResult(root=x_curr, fval=fx1, iters=1, converged=False, history=history)

    for k in range(2, max_iter + 1):
        try:
            d_prev = g_prev - x_prev
            d_curr = g_curr - x_curr

            denom = (d_curr - d_prev)
            if not is_finite(denom) or abs(denom) < 1e-14:
                w = 1.0  # iteração simples
            else:
                w = -d_curr / denom
                # clamping forte evita explosões
                w = max(w_bounds[0], min(w_bounds[1], w))

            # clamp extra
            if abs(w) > 5.0 or not math.isfinite(w):
                w = 1.0

            x_try = x_curr + w * d_curr
            if not is_finite(x_try):
                w = 1.0
                x_try = x_curr + d_curr  # iteração simples

            fx_try = f_original(x_try)

            # backtracking simples se piorar muito o residual
            if math.isfinite(fx_try) and math.isfinite(f_original(x_curr)) and abs(fx_try) > 10.0 * abs(f_original(x_curr)):
                w = 1.0
                x_try = x_curr + d_curr
                fx_try = f_original(x_try)

            # agora siga com x_new = x_try / fx_new = fx_try
            x_new, fx_new = x_try, fx_try

            # se estamos em variável y (caso B), projete para manter y>0
            if x_new <= 1e-12:
                x_new = 1e-12

            history.append((k, x_new, fx_new))

            step_ok = abs(x_new - x_curr) <= eps_abs_x
            res_ok  = abs(fx_new)        <= delta_f
            if step_ok and res_ok:
                return RootResult(root=x_new, fval=fx_new, iters=k, converged=True, history=history)

            # prepara próxima
            x_prev, x_curr = x_curr, x_new
            g_prev, g_curr = g_curr, g(x_curr)

            # se g_curr não-finito, aborta
            if not is_finite(g_curr):
                return RootResult(root=x_curr, fval=fx_new, iters=k, converged=False, history=history)

        except (ValueError, ZeroDivisionError, OverflowError):
            return RootResult(root=x_curr, fval=f_original(x_curr), iters=k-1, converged=False, history=history)

    return RootResult(root=x_curr, fval=f_original(x_curr), iters=max_iter, converged=False, history=history)


# ============================================================================
# PROBLEMA A: f(x) = exp(-x) - 2*sqrt(x) no intervalo [0,1]
# ============================================================================

def f_problema_a(x: float) -> float:
    """Função do Problema A: f(x) = exp(-x) - 2*sqrt(x)"""
    if x < 0:
        raise ValueError("x deve ser >= 0 para sqrt(x)")
    return math.exp(-x) - 2 * math.sqrt(x)


def g1_problema_a(x: float) -> float:
    """
    g1(x) = 0.25 * exp(-2*x)
    Derivada de: 2*sqrt(x) = exp(-x) => sqrt(x) = 0.5*exp(-x) => x = (0.5*exp(-x))^2
    """
    return 0.25 * math.exp(-2 * x)


def g2_problema_a(x: float) -> float:
    """
    g2(x) = -ln(2) - 0.5*ln(x)
    Derivada de: exp(-x) = 2*sqrt(x) => -x = ln(2) + 0.5*ln(x)
    """
    if x <= 0:
        raise ValueError("x deve ser > 0 para ln(x)")
    return -math.log(2) - 0.5 * math.log(x)


def solve_problema_a() -> Dict[str, RootResult]:
    """Resolve o Problema A com os três métodos."""
    print("\n" + "="*60)
    print("PROBLEMA A: f(x) = exp(-x) - 2*sqrt(x)")
    print("="*60)
    
    results = {}
    
    # Substituições Sucessivas com g1
    print("\nExecutando Substituições Sucessivas (g1)...")
    result_g1 = fixed_point(g1_problema_a, x0=0.2, f_original=f_problema_a)
    results['subst_g1'] = result_g1
    save_history_csv("problemaA_subst_g1.csv", result_g1.history, "A", "Subst_g1")
    
    # Substituições Sucessivas com g2
    print("Executando Substituições Sucessivas (g2)...")
    result_g2 = fixed_point(g2_problema_a, x0=0.5, f_original=f_problema_a)
    results['subst_g2'] = result_g2
    save_history_csv("problemaA_subst_g2.csv", result_g2.history, "A", "Subst_g2")
    
    # Wegstein com g1
    print("Executando Wegstein (g1)...")
    x1 = g1_problema_a(0.2)
    result_wegstein = wegstein(g1_problema_a, x0=0.2, x1=x1, f_original=f_problema_a)
    results['wegstein'] = result_wegstein
    save_history_csv("problemaA_wegstein.csv", result_wegstein.history, "A", "Wegstein")
    
    return results


# ============================================================================
# PROBLEMA B: Equação de escoamento turbulento
# ============================================================================

def f_problema_b(f: float, Re: float = 5000) -> float:
    """
    Função do Problema B: F(f) = 1/sqrt(f) - 1.74*ln(Re*sqrt(f)) + 0.4
    Equação de escoamento turbulento em tubo liso
    """
    if f <= 0:
        raise ValueError("f deve ser > 0")
    if Re * math.sqrt(f) <= 0:
        raise ValueError("Re*sqrt(f) deve ser > 0")
    
    return 1/math.sqrt(f) - 1.74 * math.log(Re * math.sqrt(f)) + 0.4


def g1_y_problema_b(y: float, Re: float = 5000) -> float:
    """
    g1_y(y) = 1.74*(ln(Re) - ln(y)) - 0.4
    Forma em y = 1/sqrt(f)
    """
    if y <= 0:
        raise ValueError("y deve ser > 0")
    return 1.74 * (math.log(Re) - math.log(y)) - 0.4


def g2_f_problema_b(f: float, Re: float = 5000) -> float:
    """
    g2_f(f) = [1 / (1.74*ln(Re*sqrt(f)) - 0.4)]^2
    Forma direta em f
    """
    if f <= 0:
        raise ValueError("f deve ser > 0")
    if Re * math.sqrt(f) <= 0:
        raise ValueError("Re*sqrt(f) deve ser > 0")
    
    denominator = 1.74 * math.log(Re * math.sqrt(f)) - 0.4
    if abs(denominator) < 1e-12:
        raise ValueError("Denominador muito próximo de zero")
    
    return (1 / denominator) ** 2


def f_blasius(Re: float = 5000) -> float:
    """Estimativa inicial de Blasius para f."""
    return 0.316 * (Re ** (-0.25))


def solve_problema_b() -> Dict[str, RootResult]:
    """Resolve o Problema B com os três métodos."""
    print("\n" + "="*60)
    print("PROBLEMA B: Equação de escoamento turbulento (Re=5000)")
    print("="*60)
    
    Re = 5000
    f_blas = f_blasius(Re)
    y0 = 1 / math.sqrt(f_blas)
    
    results = {}
    
    # Substituições Sucessivas com g1_y (forma em y)
    print(f"\nExecutando Substituições Sucessivas (g1_y)...")
    print(f"Chute inicial: y0 = {y0:.6f} (f_blasius = {f_blas:.6f})")
    
    def f_original_y(y: float) -> float:
        """Função original em termos de y para avaliação."""
        f = 1 / (y ** 2)
        return f_problema_b(f, Re)
    
    result_g1_y = fixed_point(g1_y_problema_b, x0=y0, f_original=f_original_y)
    # Converte resultado de y para f com verificação de overflow
    try:
        if abs(result_g1_y.root) < 1e-10:
            f_result = float('inf')
        else:
            f_result = 1 / (result_g1_y.root ** 2)
        result_g1_y.root = f_result
        result_g1_y.fval = f_problema_b(f_result, Re)
    except (OverflowError, ZeroDivisionError):
        result_g1_y.root = float('inf')
        result_g1_y.fval = float('inf')
        result_g1_y.converged = False
    results['subst_g1_y'] = result_g1_y
    save_history_csv("problemaB_subst_g1.csv", result_g1_y.history, "B", "Subst_g1_y")
    
    # Substituições Sucessivas com g2_f (forma direta em f)
    print("Executando Substituições Sucessivas (g2_f)...")
    result_g2_f = fixed_point(g2_f_problema_b, x0=f_blas, f_original=lambda f: f_problema_b(f, Re))
    results['subst_g2_f'] = result_g2_f
    save_history_csv("problemaB_subst_g2.csv", result_g2_f.history, "B", "Subst_g2_f")
    
    # Wegstein com g1_y
    print("Executando Wegstein (g1_y)...")
    y1 = g1_y_problema_b(y0)
    result_wegstein_y = wegstein(g1_y_problema_b, x0=y0, x1=y1, f_original=f_original_y)
    # Converte resultado de y para f com verificação de overflow
    try:
        if abs(result_wegstein_y.root) < 1e-10:
            f_result_weg = float('inf')
        else:
            f_result_weg = 1 / (result_wegstein_y.root ** 2)
        result_wegstein_y.root = f_result_weg
        result_wegstein_y.fval = f_problema_b(f_result_weg, Re)
    except (OverflowError, ZeroDivisionError):
        result_wegstein_y.root = float('inf')
        result_wegstein_y.fval = float('inf')
        result_wegstein_y.converged = False
    results['wegstein_y'] = result_wegstein_y
    save_history_csv("problemaB_wegstein.csv", result_wegstein_y.history, "B", "Wegstein_y")
    
    # Wegstein em f (alternativa robusta)
    print("Executando Wegstein (g2_f)...")
    f1 = g2_f_problema_b(f_blas, Re)
    result_wegstein_f = wegstein(lambda z: g2_f_problema_b(z, Re), x0=f_blas, x1=f1, f_original=lambda z: f_problema_b(z, Re))
    results['wegstein_f'] = result_wegstein_f
    save_history_csv("problemaB_wegstein_f.csv", result_wegstein_f.history, "B", "Wegstein_f")
    
    return results


def print_results_table(problem: str, results: Dict[str, RootResult]):
    """Imprime tabela de resultados para um problema."""
    print(f"\nResultados - {problem}:")
    print("-" * 70)
    print(f"{'Método':15} | {'Raiz/f':12} | {'|F(raiz)|':12} | {'Iter':4} | {'OK'}")
    print("-" * 70)
    
    for method, result in results.items():
        if method == 'wegstein_y':
            method_name = 'Wegstein Y'
        elif method == 'wegstein_f':
            method_name = 'Wegstein F'
        else:
            method_name = method.replace('_', ' ').title()
        print(format_table_row(method_name, result.root, result.fval, result.iters, result.converged))
    
    print("-" * 70)


def print_summary(results_a: Dict[str, RootResult], results_b: Dict[str, RootResult]):
    """Imprime resumo comparativo dos métodos."""
    print("\n" + "="*60)
    print("RESUMO COMPARATIVO")
    print("="*60)
    
    # Problema A
    print("\nProblema A:")
    converged_a = [k for k, v in results_a.items() if v.converged]
    if converged_a:
        fastest_a = min(converged_a, key=lambda k: results_a[k].iters)
        print(f"  - Método mais rápido: {fastest_a} ({results_a[fastest_a].iters} iterações)")
    else:
        print("  - Nenhum método convergiu")
    
    # Problema B
    print("\nProblema B:")
    converged_b = [k for k, v in results_b.items() if v.converged]
    if converged_b:
        fastest_b = min(converged_b, key=lambda k: results_b[k].iters)
        print(f"  - Método mais rápido: {fastest_b} ({results_b[fastest_b].iters} iterações)")
    else:
        print("  - Nenhum método convergiu")
    
    print("\nObservações:")
    print("- Wegstein geralmente acelera a convergência quando aplicável")
    print("- Substituições Sucessivas podem ser sensíveis ao chute inicial")
    print("- Formas diferentes de g(x) podem ter diferentes comportamentos de convergência")


if __name__ == "__main__":
    print("Lab 2 - Métodos de Raízes Não Lineares (Estudante 3)")
    print("Implementação de Substituições Sucessivas e Wegstein")
    print(f"Executado em: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Resolve Problema A
    results_a = solve_problema_a()
    print_results_table("Problema A", results_a)
    
    # Resolve Problema B
    results_b = solve_problema_b()
    print_results_table("Problema B", results_b)
    
    # Resumo comparativo
    print_summary(results_a, results_b)
    
    print(f"\nHistóricos salvos em: ./out/")
    print("Execução concluída!")
