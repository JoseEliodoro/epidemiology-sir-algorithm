
**Epidemiologia SIR — Simulação e Análise**

Projeto para simular um modelo SIR (Susceptible-Infected-Recovered) em Python.

**Resumo**:
- **Propósito:**: Simular a dinâmica de transmissão de uma doença usando o modelo SIR e gerar saídas numéricas/visuais em `data_output/`.


**Como usar**

- **Instalar dependências:**: ative seu ambiente virtual ou crie um novo e instale as dependências do `requirements.txt`.

```bash
# Instalar dependências
pip install -r requirements.txt
```

- **Executar simulação:**: execute o script principal `simulation.py`.

```bash
# Rodar a simulação
python simulation.py
```

Observação: O comportamento exato e parâmetros (taxas de transmissão/recuperação, passos de tempo, número de iterações) podem ser definidos dentro de `simulation.py` ou parâmetros de linha de comando se implementados.

**Saída**
- **Pasta de saída:**: resultados e gráficos (quando gerados) são salvos em `data_output/`.

**Estrutura do repositório**
- **`simulation.py`**: script principal que executa a simulação SIR.
- **`requirements.txt`**: dependências Python necessárias.
- **`data_output/`**: pasta onde os resultados são salvos.

