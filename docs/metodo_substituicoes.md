# Método de Substituições Sucessivas (Ponto Fixo)

## Teoria

O método de Substituições Sucessivas, também conhecido como método de ponto fixo, é uma técnica iterativa para encontrar raízes de equações não lineares.

### Algoritmo

Dada uma equação `f(x) = 0`, reescrevemos como `x = g(x)`, onde `g(x)` é uma função de iteração. O método então itera:

```
x_{k+1} = g(x_k)
```

### Critérios de Convergência

O método converge se:

1. **Passo pequeno:** `|x_{k+1} - x_k| <= eps_abs_x`
2. **Residual pequeno:** `|f(x_{k+1})| <= delta_f`

Ambos os critérios devem ser satisfeitos simultaneamente.

### Condições de Convergência

Para garantir convergência, a função `g(x)` deve satisfazer:

- `|g'(x)| < 1` em uma vizinhança da raiz
- A função deve ser contínua e diferenciável

## Implementação

### Função Principal

```python
def fixed_point(
    g: Callable[[float], float],
    x0: float,
    eps_abs_x: float = 1e-6,
    delta_f: float = 1e-6,
    max_iter: int = 500,
    f_original: Optional[Callable[[float], float]] = None
) -> RootResult:
```

### Parâmetros

- `g`: Função de iteração `g(x)`
- `x0`: Chute inicial
- `eps_abs_x`: Tolerância para o passo (padrão: 1e-6)
- `delta_f`: Tolerância para o residual (padrão: 1e-6)
- `max_iter`: Número máximo de iterações (padrão: 500)
- `f_original`: Função original `f(x)` para avaliação do residual

### Retorno

```python
@dataclass
class RootResult:
    root: float           # Raiz encontrada
    fval: float           # Valor de f(raiz)
    iters: int            # Número de iterações
    converged: bool       # Se convergiu
    history: List[Tuple[int, float, float]]  # Histórico (iter, x, f(x))
```

## Aplicações nos Problemas

### Problema A: f(x) = exp(-x) - 2\*sqrt(x)

**Funções de iteração:**

1. **g1(x) = 0.25 * exp(-2*x)**

   - Derivada de: `2*sqrt(x) = exp(-x)`
   - Chute inicial: `x0 = 0.2`
   - Comportamento: Converge rapidamente

2. **g2(x) = -ln(2) - 0.5\*ln(x)**
   - Derivada de: `exp(-x) = 2*sqrt(x)`
   - Chute inicial: `x0 = 0.5`
   - Comportamento: Pode sair do domínio (x > 0)

### Problema B: Escoamento Turbulento

**Formas de iteração:**

1. **Forma Y:** `g1_y(y) = 1.74*(ln(Re) - ln(y)) - 0.4`

   - Variável: `y = 1/sqrt(f)`
   - Chute inicial: `y0 = 1/sqrt(f_blasius)`
   - Conversão: `f = 1/y²`

2. **Forma F:** `g2_f(f) = [1/(1.74*ln(Re*sqrt(f)) - 0.4)]²`
   - Variável: `f` (direta)
   - Chute inicial: `f_blasius = 0.316 * Re^(-0.25)`

## Vantagens e Desvantagens

### Vantagens

- Simplicidade conceitual
- Fácil implementação
- Convergência rápida quando bem escolhida g(x)

### Desvantagens

- Sensível à escolha de g(x)
- Pode divergir se |g'(x)| >= 1
- Convergência linear (lenta)

## Exemplo de Uso

```python
# Definir função de iteração
def g(x):
    return 0.25 * math.exp(-2 * x)

# Definir função original
def f(x):
    return math.exp(-x) - 2 * math.sqrt(x)

# Executar método
result = fixed_point(g, x0=0.2, f_original=f)

# Verificar resultado
if result.converged:
    print(f"Raiz: {result.root}")
    print(f"f(raiz): {result.fval}")
    print(f"Iterações: {result.iters}")
```
