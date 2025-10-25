# Método de Wegstein

## Teoria

O método de Wegstein é uma técnica de aceleração aplicada sobre iterações de ponto fixo. Ele usa informações de duas iterações consecutivas para acelerar a convergência.

### Algoritmo

Dada uma função de iteração `g(x)`, o método Wegstein calcula:

1. **Diferenças:** `d_k = g(x_k) - x_k` e `d_{k-1} = g(x_{k-1}) - x_{k-1}`
2. **Fator de aceleração:** `w = -d_k / (d_k - d_{k-1})`
3. **Nova iteração:** `x_{k+1} = x_k + w * d_k`

### Interpretação

- Se `w = 1`: Volta à iteração simples `x_{k+1} = g(x_k)`
- Se `w = 0`: Mantém `x_{k+1} = x_k` (sem progresso)
- Se `w > 1`: Acelera na direção do movimento
- Se `w < 0`: Inverte a direção

## Implementação Robusta

### Salvaguardas Implementadas

1. **Clamp do fator w:**

   ```python
   w = max(w_bounds[0], min(w_bounds[1], w))  # Padrão: [-10, 10]
   ```

2. **Clamp extra:**

   ```python
   if abs(w) > 5.0 or not math.isfinite(w):
       w = 1.0  # Volta à iteração simples
   ```

3. **Backtracking:**

   ```python
   if abs(fx_try) > 10.0 * abs(f_original(x_curr)):
       w = 1.0  # Re-tenta com iteração simples
   ```

4. **Projeção de domínio:**
   ```python
   if x_new <= 1e-12:  # Para variáveis que devem ser > 0
       x_new = 1e-12
   ```

### Função Principal

```python
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
```

### Parâmetros

- `g`: Função de iteração `g(x)`
- `x0, x1`: Dois chutes iniciais
- `eps_abs_x`: Tolerância para o passo
- `delta_f`: Tolerância para o residual
- `max_iter`: Número máximo de iterações
- `w_bounds`: Limites para o fator w (padrão: [-10, 10])
- `f_original`: Função original `f(x)` para avaliação

## Aplicações nos Problemas

### Problema A: f(x) = exp(-x) - 2\*sqrt(x)

**Configuração:**

- Função base: `g1(x) = 0.25 * exp(-2*x)`
- Chutes: `x0 = 0.2`, `x1 = g1(x0)`
- Comportamento: Pode não convergir devido à natureza da função

### Problema B: Escoamento Turbulento

**Duas implementações:**

1. **Wegstein Y:**

   - Base: `g1_y(y) = 1.74*(ln(Re) - ln(y)) - 0.4`
   - Variável: `y = 1/sqrt(f)`
   - Conversão: `f = 1/y²`

2. **Wegstein F:**
   - Base: `g2_f(f) = [1/(1.74*ln(Re*sqrt(f)) - 0.4)]²`
   - Variável: `f` (direta)
   - Mais robusta numericamente

## Vantagens e Desvantagens

### Vantagens

- Acelera convergência quando aplicável
- Fallback automático para iteração simples
- Robustez numérica com salvaguardas
- Funciona bem com funções bem comportadas

### Desvantagens

- Mais complexo que ponto fixo simples
- Pode ser instável com funções mal condicionadas
- Requer dois chutes iniciais
- Overhead computacional adicional

## Exemplo de Uso

```python
# Definir função de iteração
def g(x):
    return 1.74 * (math.log(5000) - math.log(x)) - 0.4

# Definir função original
def f(x):
    y = x
    f_val = 1/(y**2)
    return 1/math.sqrt(f_val) - 1.74 * math.log(5000 * math.sqrt(f_val)) + 0.4

# Executar Wegstein
x1 = g(x0)
result = wegstein(g, x0=x0, x1=x1, f_original=f)

# Verificar resultado
if result.converged:
    print(f"Raiz: {result.root}")
    print(f"f(raiz): {result.fval}")
    print(f"Iterações: {result.iters}")
```

## Comparação com Ponto Fixo

| Aspecto      | Ponto Fixo      | Wegstein       |
| ------------ | --------------- | -------------- |
| Simplicidade | Alta            | Média          |
| Velocidade   | Linear          | Acelerada      |
| Robustez     | Baixa           | Alta           |
| Overhead     | Baixo           | Médio          |
| Convergência | Sensível a g(x) | Mais tolerante |

## Dicas de Uso

1. **Escolha bons chutes iniciais:** Próximos da raiz esperada
2. **Monitore o fator w:** Valores extremos indicam problemas
3. **Use salvaguardas:** Sempre ative clamp e backtracking
4. **Teste diferentes formas:** Y vs F podem ter comportamentos diferentes
5. **Verifique domínio:** Garanta que variáveis permanecem válidas
