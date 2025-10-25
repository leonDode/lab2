# Problemas Teóricos - Métodos de Raízes Não Lineares

## Problema A: Equação Transcendental

### Enunciado

Encontrar a raiz positiva da equação:

```
f(x) = exp(-x) - 2*sqrt(x) = 0
```

no intervalo [0, 1].

### Análise Matemática

**Domínio:** x ≥ 0 (devido ao sqrt(x))

**Comportamento:**

- f(0) = exp(0) - 2\*sqrt(0) = 1 - 0 = 1
- f(1) = exp(-1) - 2\*sqrt(1) ≈ 0.368 - 2 = -1.632
- Pelo teorema do valor intermediário, existe raiz em [0,1]

**Derivada:**

```
f'(x) = -exp(-x) - 1/sqrt(x)
```

A derivada é sempre negativa, indicando função decrescente.

### Funções de Iteração

#### g1(x) = 0.25 * exp(-2*x)

**Derivação:**

```
2*sqrt(x) = exp(-x)
sqrt(x) = 0.5*exp(-x)
x = (0.5*exp(-x))² = 0.25*exp(-2*x)
```

**Análise de convergência:**

```
g1'(x) = 0.25*(-2)*exp(-2*x) = -0.5*exp(-2*x)
|g1'(x)| = 0.5*exp(-2*x) < 0.5 < 1
```

Convergência garantida.

#### g2(x) = -ln(2) - 0.5\*ln(x)

**Derivação:**

```
exp(-x) = 2*sqrt(x)
-x = ln(2) + 0.5*ln(x)
x = -ln(2) - 0.5*ln(x)
```

**Análise de convergência:**

```
g2'(x) = -0.5/x
|g2'(x)| = 0.5/|x|
```

Para x próximo de 0, |g2'(x)| pode ser > 1, causando divergência.

### Resultados Esperados

- **g1:** Converge para x ≈ 0.175867
- **g2:** Pode divergir ou sair do domínio (x > 0)

## Problema B: Escoamento Turbulento

### Enunciado

Resolver a equação implícita para o fator de atrito f em escoamento turbulento em tubo liso:

```
1/√f = 1.74*ln(Re*√f) - 0.4
```

onde Re = 5000 (número de Reynolds).

### Contexto Físico

**Equação de Colebrook-White:**

```
1/√f = -2.0*log10(2.51/(Re*√f) + k/(3.7*D))
```

Para tubo liso (k = 0):

```
1/√f = -2.0*log10(2.51/(Re*√f))
```

Convertendo para ln e ajustando constantes:

```
1/√f = 1.74*ln(Re*√f) - 0.4
```

### Estimativa Inicial (Blasius)

```
f_blasius = 0.316 * Re^(-0.25)
f_blasius = 0.316 * 5000^(-0.25) ≈ 0.0376
```

### Formas de Iteração

#### Forma Y: g1_y(y) = 1.74\*(ln(Re) - ln(y)) - 0.4

**Variável:** y = 1/√f

**Derivação:**

```
1/√f = 1.74*ln(Re*√f) - 0.4
y = 1.74*ln(Re/y) - 0.4
y = 1.74*(ln(Re) - ln(y)) - 0.4
```

**Chute inicial:** y₀ = 1/√f_blasius ≈ 5.16

#### Forma F: g2_f(f) = [1/(1.74*ln(Re*√f) - 0.4)]²

**Variável:** f (direta)

**Derivação:**

```
1/√f = 1.74*ln(Re*√f) - 0.4
√f = 1/(1.74*ln(Re*√f) - 0.4)
f = [1/(1.74*ln(Re*√f) - 0.4)]²
```

**Chute inicial:** f₀ = f_blasius ≈ 0.0376

### Análise Numérica

**Domínio:** f > 0, Re\*√f > 0

**Comportamento esperado:**

- f ≈ 0.00933 (valor típico para Re = 5000)
- Ambas as formas devem convergir
- Forma F é mais robusta numericamente

### Resultados Esperados

- **g1_y:** Converge para y ≈ 10.35, f ≈ 0.00933
- **g2_f:** Converge para f ≈ 0.00933
- **Wegstein:** Acelera convergência em ambas as formas

## Comparação dos Problemas

| Aspecto          | Problema A      | Problema B        |
| ---------------- | --------------- | ----------------- |
| **Tipo**         | Transcendental  | Físico/Engenharia |
| **Domínio**      | [0,1]           | f > 0             |
| **Complexidade** | Baixa           | Média             |
| **Convergência** | Sensível a g(x) | Mais robusta      |
| **Aplicação**    | Matemática pura | Engenharia        |

## Dicas para Implementação

### Problema A

1. **g1 é mais confiável** que g2
2. **Monitore domínio** para g2 (x > 0)
3. **Chutes iniciais** próximos à raiz esperada

### Problema B

1. **Forma F é mais robusta** que forma Y
2. **Use estimativa de Blasius** como chute inicial
3. **Monitore overflow** na conversão y ↔ f
4. **Wegstein funciona bem** em ambas as formas

## Validação dos Resultados

### Problema A

- Verificar: f(raiz) ≈ 0
- Raiz esperada: x ≈ 0.175867

### Problema B

- Verificar: F(f) ≈ 0
- f esperado: ≈ 0.00933
- Validar fisicamente: f deve estar em [0.01, 0.1] para Re = 5000
