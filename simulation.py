import numpy as np
import matplotlib.pyplot as plt
import os

# --- Parâmetros de Simulação ---
# N: População Total
# I0: Infectados Iniciais
# R0_param: Recuperados Iniciais
# beta (β): Taxa de Transmissão (probabilidade de infecção por contato)
# gamma (γ): Taxa de Recuperação (1 / duração média da infecção)
# total_time: Tempo total de simulação (dias)
# time_steps: Número de passos de tempo (igual a total_time para dt=1)

# Cenários de Doenças para Comparação
datas = [
    {
        "nome": "Cenário 1: Doença Altamente Contagiosa (Alta beta)",
        "N": 100000,
        "I0": 10,
        "R0_param": 0,
        "beta": 0.5, # Taxa de transmissão alta
        "gamma": 0.1, # Taxa de recuperação moderada (R0 = 5.0)
        "total_time": 200,
        "time_steps": 200,
        "cor": "red",
    },
    {
        "nome": "Cenário 2: Doença Levemente Contagiosa (Baixa beta)",
        "N": 100000,
        "I0": 10,
        "R0_param": 0,
        "beta": 0.2, # Taxa de transmissão baixa
        "gamma": 0.1, # Taxa de recuperação moderada (R0 = 2.0)
        "total_time": 200,
        "time_steps": 200,
        "cor": "blue",
    },
    {
        "nome": "Cenário 3: Intervenção Tardia (Redução de beta)",
        "N": 100000,
        "I0": 10,
        "R0_param": 0,
        "beta": 0.5, # Começa alta
        "gamma": 0.1,
        "intervencao_dia": 50, # Intervenção começa no dia 50
        "beta_pos_intervencao": 0.15, # Cai drasticamente
        "total_time": 200,
        "time_steps": 200,
        "cor": "green",
    },
    {
        "nome": "Cenário 4: Recuperação Lenta (Baixa gamma)",
        "N": 100000,
        "I0": 10,
        "R0_param": 0,
        "beta": 0.3, 
        "gamma": 0.05, # Taxa de recuperação muito lenta
        "total_time": 200,
        "time_steps": 200,
        "cor": "purple",
    }
]

def run_sir_model(N, I0, R0_param, beta, gamma, total_time, time_steps, intervencao_dia=None, beta_pos_intervencao=None):
    """
    Executa a simulação do Modelo SIR usando o método de Euler.

    Argumentos:
        N: População total.
        I0: Número inicial de Infectados.
        R0_param: Número inicial de Recuperados.
        beta: Taxa de transmissão (β).
        gamma: Taxa de recuperação (γ).
        total_time: Duração da simulação.
        time_steps: Número de passos de tempo.
        intervencao_dia: (Opcional) Dia em que beta muda.
        beta_pos_intervencao: (Opcional) Nova taxa beta após a intervenção.

    Retorna:
        t: Array de tempo.
        S: Array de Suscetíveis.
        I: Array de Infectados.
        R: Array de Recuperados.
    """
    
    # 1. Inicialização
    dt = total_time / time_steps
    t = np.linspace(0, total_time, time_steps + 1)
    
    # Condições iniciais
    I = np.zeros(time_steps + 1)
    R = np.zeros(time_steps + 1)
    S = np.zeros(time_steps + 1)

    I[0] = I0
    R[0] = R0_param
    S[0] = N - I0 - R0_param

    # 2. Resolução do sistema de equações (Euler Forward)
    for i in range(time_steps):
        # Verifica a intervenção para mudar beta dinamicamente
        current_beta = beta
        if intervencao_dia is not None and t[i] >= intervencao_dia:
            current_beta = beta_pos_intervencao
            
        # Cálculo das mudanças (Derivadas)
        # dS/dt = - beta * S * I / N
        d_S = -current_beta * S[i] * I[i] / N
        # dI/dt = (beta * S * I / N) - gamma * I
        d_I = (current_beta * S[i] * I[i] / N) - gamma * I[i]
        # dR/dt = gamma * I
        d_R = gamma * I[i]

        # Atualização dos compartimentos para o próximo passo
        S[i+1] = S[i] + d_S * dt
        I[i+1] = I[i] + d_I * dt
        R[i+1] = R[i] + d_R * dt
        
        # Garantir que as populações não sejam negativas (correção numérica)
        S[i+1] = max(0, S[i+1])
        I[i+1] = max(0, I[i+1])
        R[i+1] = min(N, R[i+1]) # R não pode ser maior que N

        # Normalização (S+I+R = N) devido a erros numéricos no Euler
        total = S[i+1] + I[i+1] + R[i+1]
        if total != N and total > 0:
             fator = N / total
             S[i+1] *= fator
             I[i+1] *= fator
             R[i+1] *= fator


    return t, S, I, R

def plot_single_sir_curve(t, S, I, R, data, output_dir, filename):
    """
    Plota as curvas S, I e R para um único cenário (gráfico individual).
    """
    plt.figure(figsize=(10, 6))
    
    # Calcular R0 inicial
    R0_calc = data["beta"] / data["gamma"]
    
    plt.plot(t, S, color='blue', linewidth=2, label='Suscetíveis (S)')
    plt.plot(t, I, color='red', linewidth=2, label='Infectados (I)')
    plt.plot(t, R, color='green', linewidth=2, label='Recuperados (R)')
    
    # Adicionar linha vertical para intervenção, se aplicável
    if data.get("intervencao_dia"):
        plt.axvline(x=data["intervencao_dia"], color='orange', linestyle='--', 
                    label='Início da Intervenção')

    plt.title(f'Modelo SIR: {data["nome"]} (R0 inicial ≈ {R0_calc:.1f})', fontsize=14)
    plt.xlabel('Tempo (Dias)', fontsize=12)
    plt.ylabel('População', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='center right')
    plt.ylim(bottom=0, top=data["N"] * 1.05)
    plt.xlim(left=0)
    
    # Salvar a figura
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath)
    plt.close()
    print(f"Gráfico individual salvo em: {filepath}")

def plot_sir_comparison(results, N, output_dir, filename="sir_comparison.png"):
    """
    Plota as curvas de Infectados (I) de múltiplos cenários para comparação.
    """
    plt.figure(figsize=(12, 6))
    
    # Plotar apenas a curva I (Infectados) para comparação
    for result in results:
        t, S, I, R = result["t"], result["S"], result["I"], result["R"]
        label = result["nome"]
        cor = result["cor"]
        
        # Calcular o R0 inicial para o título
        beta_inicial = result["beta"]
        gamma = result["gamma"]
        R0_calc = beta_inicial / gamma
        
        plt.plot(t, I, color=cor, linewidth=3, 
                 label=f'{label} ($R_0$ inicial ≈ {R0_calc:.1f})')
        
        # Marcar o pico de infecção
        pico_I = np.max(I)
        dia_pico = t[np.argmax(I)]
        plt.plot(dia_pico, pico_I, marker='o', color=cor, markersize=5)
        plt.annotate(f'Pico: {pico_I:,.0f}', (dia_pico, pico_I * 1.05), color=cor, 
                     fontsize=9, ha='center')

    plt.title(f'Comparação de Cenários de Contágio (População Total N={N:,})', fontsize=14)
    plt.xlabel('Tempo (Dias)', fontsize=12)
    plt.ylabel('Número de Indivíduos Infectados (I)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='upper right')
    plt.ylim(bottom=0)
    plt.xlim(left=0)
    
    # Salvar a figura
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath)
    plt.close() # Fechar a figura para liberar memória
    print(f"Gráfico de comparação salvo em: {filepath}")


# --- Execução Principal ---
if __name__ == '__main__':
    # Diretório de saída
    output_dir = "data_output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    all_results = []
    
    print("Iniciando simulações para os diferentes cenários...")

    for i, data in enumerate(datas):
        print(f"Simulando: {data['nome']}")
        
        # Verifica se há parâmetros de intervenção
        intervencao_dia = data.get("intervencao_dia")
        beta_pos_intervencao = data.get("beta_pos_intervencao")
        
        t, S, I, R = run_sir_model(
            data["N"], 
            data["I0"], 
            data["R0_param"], 
            data["beta"], 
            data["gamma"], 
            data["total_time"], 
            data["time_steps"],
            intervencao_dia,
            beta_pos_intervencao
        )
        
        # Armazena os resultados junto com metadados para o plot de comparação
        all_results.append({
            "t": t, "S": S, "I": I, "R": R, 
            "nome": data["nome"], 
            "cor": data["cor"], 
            "beta": data["beta"],
            "gamma": data["gamma"]
        })
        
        # NOVO: Plotar o gráfico individual S-I-R para este cenário
        plot_single_sir_curve(t, S, I, R, data, output_dir, f"sir_cenario_{i+1}_individual.png")

    # Plotar todos os cenários de Infectados em um único gráfico (comparação)
    N_total = datas[0]["N"]
    plot_sir_comparison(all_results, N_total, output_dir, "comparacao_sir_final.png")

    print("\nSimulações concluídas. Foram gerados 4 gráficos individuais (S-I-R) e 1 gráfico de comparação (Curva I) no diretório 'data_output'.")