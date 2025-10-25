#!/usr/bin/env python3
"""
Lab 2 - Problema A: Equação Transcendental
Estudante 3: Implementação do método de Substituições Sucessivas

Este módulo resolve o Problema A: f(x) = exp(-x) - 2*sqrt(x) = 0
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
    """Resolve o Problema A com Substituições Sucessivas."""
    print("\n" + "="*60)
    print("PROBLEMA A: f(x) = exp(-x) - 2*sqrt(x)")
    print("Método: Substituições Sucessivas (Ponto Fixo)")
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
    
    return results


def print_results_table(results: Dict[str, RootResult]):
    """Imprime tabela de resultados para o Problema A."""
    print(f"\nResultados - Problema A:")
    print("-" * 70)
    print(f"{'Método':15} | {'Raiz':12} | {'|f(raiz)|':12} | {'Iter':4} | {'OK'}")
    print("-" * 70)
    
    for method, result in results.items():
        method_name = method.replace('_', ' ').title()
        print(format_table_row(method_name, result.root, result.fval, result.iters, result.converged))
    
    print("-" * 70)


def print_summary(results: Dict[str, RootResult]):
    """Imprime resumo comparativo do Problema A."""
    print("\n" + "="*60)
    print("RESUMO COMPARATIVO - PROBLEMA A")
    print("="*60)
    
    converged_methods = [k for k, v in results.items() if v.converged]
    if converged_methods:
        fastest = min(converged_methods, key=lambda k: results[k].iters)
        print(f"\nMétodo mais rápido: {fastest} ({results[fastest].iters} iterações)")
    else:
        print("\nNenhum método convergiu")
    
    print("\nAnálise dos métodos:")
    print("- g1(x) = 0.25*exp(-2*x): Converge rapidamente (derivada < 1)")
    print("- g2(x) = -ln(2) - 0.5*ln(x): Pode divergir (derivada pode ser > 1)")
    
    print("\nObservações:")
    print("- Substituições Sucessivas são sensíveis à escolha de g(x)")
    print("- g1 é mais confiável que g2 para este problema")
    print("- Raiz esperada: x ≈ 0.175867")


if __name__ == "__main__":
    print("Lab 2 - Problema A: Equação Transcendental")
    print("Método: Substituições Sucessivas (Ponto Fixo)")
    print(f"Executado em: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Resolve Problema A
    results = solve_problema_a()
    print_results_table(results)
    
    # Resumo comparativo
    print_summary(results)
    
    print(f"\nHistóricos salvos em: ./out/")
    print("Execução concluída!")
