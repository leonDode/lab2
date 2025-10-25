#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lab 2 — Método de Wegstein (APENAS PROBLEMA B)
Estudante 3: Implementação do método de Wegstein com aceleração

Este módulo resolve o PROBLEMA B (escoamento turbulento em tubo liso)
usando o método de Wegstein, mantendo o padrão de saída dos demais arquivos.

Observações:
- Usamos ln natural.
- O f calculado é o fator de atrito de FANNING; para converter a DARCY, usar f_D = 4*f_F.
"""

import math
import time
import csv
from dataclasses import dataclass
from typing import List, Tuple, Dict, Callable, Optional
from pathlib import Path


# =============================================================================
# Infraestrutura comum
# =============================================================================

@dataclass
class RootResult:
    """Resultado de um método de busca de raiz / ponto fixo."""
    root: float            # aqui será f (Fanning) ao final
    fval: float            # |F(f)| (resíduo na raiz reportada)
    iters: int
    converged: bool
    history: List[Tuple[int, float, float, float]]  # (iter, estado_xk, delta_x, F(f_k))


def is_finite(x: float) -> bool:
    return math.isfinite(x)


def _fmt_num(x: float) -> str:
    if x is None or not math.isfinite(x):
        return f"{x}"
    return f"{x:12.6f}"


def _fmt_res(r: float) -> str:
    if r is None or not math.isfinite(r):
        return f"{r}"
    return f"{abs(r):12.2e}"


def format_table_row(method: str, root: float, fval: float, iters: int, converged: bool) -> str:
    """Linhas da tabela: se não convergir, mostra '—' em f e |F|."""
    status = "✓" if converged else "✗"
    root_s = _fmt_num(root) if converged else f"{'—':>12}"
    fval_s = _fmt_res(fval) if converged else f"{'—':>12}"
    return f"{method:15} | {root_s} | {fval_s} | {iters:4d} | {status}"


def save_history_csv(filename: str, history: List[Tuple[int, float, float, float]], problem: str, method: str, state_label: str):
    """Salva histórico de iterações em CSV."""
    Path("./out").mkdir(exist_ok=True)
    filepath = f"./out/{filename}"
    with open(filepath, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Iteracao', state_label, 'Delta_'+state_label, 'F(f_k)', 'Problema', 'Metodo'])
        for iter_num, x, dx, fx in history:
            writer.writerow([iter_num, x, dx, fx, problem, method])


# =============================================================================
# Wegstein (aceleração de ponto fixo)
# =============================================================================

def wegstein(
    g: Callable[[float], float],
    x0: float,
    x1: float,
    eps_abs_x: float = 1e-6,
    delta_f: float = 1e-6,
    max_iter: int = 500,
    w_bounds: Tuple[float, float] = (0.0, 2.0),   # clipping conservador
    accel_start_iter: int = 3,                    # só acelera a partir da 3ª iteração
    f_residual: Optional[Callable[[float], float]] = None,
    to_f: Optional[Callable[[float], float]] = None,
    domain_ok: Optional[Callable[[float], bool]] = None,
    state_label: str = "x"
) -> RootResult:
    """
    Aceleração de Wegstein sobre a iteração de ponto fixo x_{k+1} = g(x_k).

    Implementa salvaguardas recomendadas:
      - Atraso para iniciar aceleração (accel_start_iter)
      - Clipping de w em [0, 2] (ou como desejado)
      - Fallback para passo "plain" se |F| piorar
      - Projeção leve para manter domínio (via domain_ok)
      - Critério de parada combinado: |Δx| ≤ eps_abs_x e |F| ≤ delta_f (enunciado)
    """
    assert f_residual is not None, "f_residual (F(f)) é obrigatório"
    assert to_f is not None, "to_f (estado -> f) é obrigatório"

    history: List[Tuple[int, float, float, float]] = []

    # Estado inicial
    if not is_finite(x1):
        x1 = g(x0)

    # Iteração 0 (estado x0)
    f0 = f_residual(to_f(x0))
    history.append((0, x0, float('nan'), f0))

    # Iteração 1 (estado x1, sem aceleração)
    f1 = f_residual(to_f(x1))
    history.append((1, x1, abs(x1 - x0), f1))

    x_prev, x_curr = x0, x1
    try:
        g_prev = g(x_prev)
        g_curr = g(x_curr)
    except Exception:
        return RootResult(root=to_f(x_curr), fval=f_residual(to_f(x_curr)), iters=1, converged=False, history=history)

    for k in range(2, max_iter + 1):
        # Próximo valor "cru"
        x_next_plain = g(x_curr)

        # Aceleração (a partir de accel_start_iter)
        use_accel = (k >= accel_start_iter)
        if use_accel:
            d_prev = g_prev - x_prev
            d_curr = g_curr - x_curr
            denom = (d_curr - d_prev)
            if (not is_finite(denom)) or abs(denom) < 1e-14:
                w = 1.0
            else:
                # Fórmula de Wegstein tradicional
                w = -d_curr / denom
                w = max(w_bounds[0], min(w_bounds[1], w))
        else:
            w = 1.0

        # Candidato acelerado
        x_try = x_curr + w * (g_curr - x_curr)
        if (domain_ok is not None) and (not domain_ok(x_try)):
            # sai do domínio -> usa passo plain
            x_try = x_next_plain
            w = 1.0

        # Métricas
        f_curr = f_residual(to_f(x_curr))
        f_try  = f_residual(to_f(x_try))

        # Fallback se resíduo piorar muito
        if (not math.isnan(f_try)) and (not math.isinf(f_try)) and abs(f_try) > 10.0 * abs(f_curr):
            x_try = x_next_plain
            w = 1.0
            f_try = f_residual(to_f(x_try))

        dx = abs(x_try - x_curr)
        history.append((k, x_try, dx, f_try))

        # Critério de parada combinado
        if dx <= eps_abs_x and abs(f_try) <= delta_f:
            return RootResult(root=to_f(x_try), fval=f_try, iters=k, converged=True, history=history)

        # Avançar
        x_prev, x_curr = x_curr, x_try
        g_prev, g_curr = g_curr, g(x_curr)

        # Se g_curr explodir, retorna melhor estimativa atual
        if not is_finite(g_curr):
            return RootResult(root=to_f(x_curr), fval=f_residual(to_f(x_curr)), iters=k, converged=False, history=history)

    # k_max atingido
    return RootResult(root=to_f(x_curr), fval=f_residual(to_f(x_curr)), iters=max_iter, converged=False, history=history)


# =============================================================================
# PROBLEMA B: Equação de escoamento turbulento (Re=5000)
#   F(f) = 1/sqrt(f) - 1.74*ln(Re*sqrt(f)) + 0.4 = 0
#   Duas formas de iteração:
#     (1) Em y = 1/sqrt(f): y = 1.74*(ln Re - ln y) - 0.4
#     (2) Em f diretamente: f = [1 / (1.74*ln(Re*sqrt(f)) - 0.4)]^2
# =============================================================================

def F_problema_b(f: float, Re: float = 5000) -> float:
    if f <= 0:
        raise ValueError("f deve ser > 0")
    val = Re * math.sqrt(f)
    if val <= 0:
        raise ValueError("Re*sqrt(f) deve ser > 0")
    return 1.0 / math.sqrt(f) - 1.74 * math.log(val) + 0.4  # ln natural


def g1_y_problema_b(y: float, Re: float = 5000) -> float:
    if y <= 0:
        raise ValueError("y deve ser > 0")
    return 1.74 * (math.log(Re) - math.log(y)) - 0.4   # y_{k+1} = 1.74*ln(Re/y_k) - 0.4


def g2_f_problema_b(f: float, Re: float = 5000) -> float:
    if f <= 0:
        raise ValueError("f deve ser > 0")
    den = 1.74 * math.log(Re * math.sqrt(f)) - 0.4
    if abs(den) < 1e-12:
        raise ValueError("Denominador ~ 0 em g2_f")
    return (1.0 / den) ** 2


def f_blasius(Re: float = 5000) -> float:
    """Estimativa inicial de Blasius (Darcy): f_D ≈ 0.316 * Re^{-0.25}."""
    return 0.316 * (Re ** (-0.25))


def solve_problema_b() -> Dict[str, RootResult]:
    """Resolve o Problema B com Wegstein (duas formas)."""
    print("Lab 2 - Método de Wegstein")
    print("Estudante 3: Implementação do método de Wegstein com aceleração")
    print(f"Executado em: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    print("="*60)
    print("PROBLEMA B: Equação de escoamento turbulento (Re=5000)")
    print("Método: Wegstein com Aceleração")
    print(f"Executado em: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*60)

    # Configuração
    Re = 5000
    eps_abs = 1e-6
    delta    = 1e-6
    k_max    = 500

    # Chute de Blasius (DARCy) e conversão para y
    f0_blas_darcy = f_blasius(Re)
    y0            = 1.0 / math.sqrt(f0_blas_darcy)

    results: Dict[str, RootResult] = {}

    # ---------------------------
    # Wegstein em y (forma g1_y)
    # ---------------------------
    print("\nExecutando Wegstein (g1_y)...")
    print(f"Chute inicial: y0 = {y0:.6f} (f_blasius_darcy = {f0_blas_darcy:.6f})")

    res_y = wegstein(
        g=lambda t: g1_y_problema_b(t, Re),
        x0=y0,
        x1=g1_y_problema_b(y0, Re),
        eps_abs_x=eps_abs,
        delta_f=delta,
        max_iter=k_max,
        w_bounds=(0.0, 2.0),
        accel_start_iter=3,
        f_residual=lambda y: F_problema_b(1.0/(y*y), Re),
        to_f=lambda y: 1.0/(y*y),
        domain_ok=lambda y: y > 0.0,
        state_label="y"
    )
    results['wegstein_y'] = res_y
    save_history_csv("problemaB_wegstein_y.csv", res_y.history, "B", "Wegstein_y", state_label="y")

    # ---------------------------
    # Wegstein em f (forma g2_f)
    # ---------------------------
    print("\nExecutando Wegstein (g2_f)...")
    print(f"Chute inicial: f0 = {f0_blas_darcy:.6f} (mesmo valor de Blasius/Darcy)")

    res_f = wegstein(
        g=lambda z: g2_f_problema_b(z, Re),
        x0=f0_blas_darcy,
        x1=g2_f_problema_b(f0_blas_darcy, Re),
        eps_abs_x=eps_abs,
        delta_f=delta,
        max_iter=k_max,
        w_bounds=(0.0, 2.0),
        accel_start_iter=3,
        f_residual=lambda f: F_problema_b(f, Re),
        to_f=lambda f: f,
        domain_ok=lambda f: f > 0.0,
        state_label="f"
    )
    results['wegstein_f'] = res_f
    save_history_csv("problemaB_wegstein_f.csv", res_f.history, "B", "Wegstein_f", state_label="f")

    # ---------------------------
    # Impressão dos resultados
    # ---------------------------
    print_results_table("Problema B", results)

    # Conversão para Darcy (contexto físico)
    print("\nConversão (contexto físico): ln natural; f na tabela é FANNING.")
    def to_darcy(f_fan: float) -> float: return 4.0 * f_fan

    # Escolhe um f convergente para exibir conversão
    if results['wegstein_y'].converged:
        f_fan = results['wegstein_y'].root
    elif results['wegstein_f'].converged:
        f_fan = results['wegstein_f'].root
    else:
        f_fan = float('nan')

    if math.isfinite(f_fan):
        print(f"f_Fanning ≈ {f_fan:.6f}  ->  f_Darcy = 4*f_Fanning ≈ {to_darcy(f_fan):.6f}")

    print(f"\nHistóricos salvos em: ./out/")
    print("Execução concluída!\n")

    return results


def print_results_table(problem: str, results: Dict[str, RootResult]):
    print(f"\nResultados - {problem}:")
    print("-" * 70)
    print(f"{'Método':15} | {'f (raiz)':12} | {'|F(f)|':12} | {'Iter':4} | {'OK'}")
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


# =============================================================================
# Execução direta (APENAS Problema B)
# =============================================================================
if __name__ == "__main__":
    solve_problema_b()
