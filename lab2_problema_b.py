#!/usr/bin/env python3
"""
Lab 2 - Problema B: Escoamento Turbulento
Estudante 3: Implementação do método de Substituições Sucessivas

Este módulo resolve o Problema B: Equação de escoamento turbulento em tubo liso
usando o método de Substituições Sucessivas (Ponto Fixo).
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
    """Resolve o Problema B com Substituições Sucessivas."""
    print("\n" + "="*60)
    print("PROBLEMA B: Equação de escoamento turbulento (Re=5000)")
    print("Método: Substituições Sucessivas (Ponto Fixo)")
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
    
    return results


def print_results_table(results: Dict[str, RootResult]):
    """Imprime tabela de resultados para o Problema B."""
    print(f"\nResultados - Problema B:")
    print("-" * 70)
    print(f"{'Método':15} | {'f':12} | {'|F(f)|':12} | {'Iter':4} | {'OK'}")
    print("-" * 70)
    
    for method, result in results.items():
        method_name = method.replace('_', ' ').title()
        print(format_table_row(method_name, result.root, result.fval, result.iters, result.converged))
    
    print("-" * 70)


def print_summary(results: Dict[str, RootResult]):
    """Imprime resumo comparativo do Problema B."""
    print("\n" + "="*60)
    print("RESUMO COMPARATIVO - PROBLEMA B")
    print("="*60)
    
    converged_methods = [k for k, v in results.items() if v.converged]
    if converged_methods:
        fastest = min(converged_methods, key=lambda k: results[k].iters)
        print(f"\nMétodo mais rápido: {fastest} ({results[fastest].iters} iterações)")
    else:
        print("\nNenhum método convergiu")
    
    print("\nAnálise dos métodos:")
    print("- g1_y(y): Forma em y = 1/sqrt(f), convergência estável")
    print("- g2_f(f): Forma direta em f, mais robusta numericamente")
    
    print("\nContexto físico:")
    print("- Re = 5000 (número de Reynolds)")
    print("- f_blasius ≈ 0.0376 (estimativa inicial)")
    print("- f esperado ≈ 0.00933 (fator de atrito)")
    
    print("\nObservações:")
    print("- Ambas as formas convergem para o mesmo resultado")
    print("- Forma F é mais robusta para implementação numérica")
    print("- Resultado fisicamente válido para escoamento turbulento")


if __name__ == "__main__":
    print("Lab 2 - Problema B: Escoamento Turbulento")
    print("Método: Substituições Sucessivas (Ponto Fixo)")
    print(f"Executado em: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Resolve Problema B
    results = solve_problema_b()
    print_results_table(results)
    
    # Resumo comparativo
    print_summary(results)
    
    print(f"\nHistóricos salvos em: ./out/")
    print("Execução concluída!")
