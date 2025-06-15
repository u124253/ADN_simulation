import random, math, numpy as np, matplotlib.pyplot as plt
from datetime import datetime
def Pe_teorico(p):
    """Probabilidad de error teórica MV / ML (n = 3)."""
    if p<= 0.75:
        return (7*p**2 - 4*p**3)/3
    
    else:
        return 1 - 16 * (p / 3)**3

def N_requerido(p, k=200):
    """Devuelve N = k / Pe, redondeado a 10^m más cercano."""
    Pe = Pe_teorico(p)
    N0  = k / Pe
    
    # redondear a potencia de 10
    #pot10 = 10 ** round(math.log10(N))
    #Redondear al multiplo mas proximo de 256
    
    
    return max(1, int(N0 ))

def log_ml_decoder(recibido, p):
    p = max(p, 1e-15)
    scores_log = {}
    max_log = float("-inf")

    for nuc in prior.keys():
        conteo = recibido.count(nuc)
        # Protección contra log(0)
        if p == 0:
            prob_log = float('-inf')  # log(0) es -infinito
        else:
            prob_log = conteo * math.log(1 - p) + (n - conteo) * math.log(p / 3)

        scores_log[nuc] = prob_log
        if prob_log > max_log:
            max_log = prob_log

    epsilon = 1e-100
    max_nucs = [k for k, v in scores_log.items() if abs(v - max_log) < epsilon]
    return max_nucs[0] if len(max_nucs) == 1 else random.choice(max_nucs)

    # Obtener los nucleótidos con máxima log-probabilidad
 # Comparación de máximos con tolerancia
    epsilon = 1e-100
    max_nucs = [k for k, v in scores_log.items() if abs(v - max_log) < epsilon]
    # Si hay un único ganador, lo devolvemos; si no, elegimos aleatoriamente entre los mejores
    return max_nucs[0] if len(max_nucs) == 1 else random.choice(max_nucs)

def rand_nucleotido(prior):
    nucleotidos = list(prior.keys())
    pesos = list(prior.values())
    return random.choices(nucleotidos, weights=pesos, k=1)[0]

def rand_prob():
    return random.random()

"""
En caso de que ocurra el cambio de letra se selecciona alguna de las 3 restantes diferentes a la original
"""
def rand_nuc_diferente(nucleotido):
    #Selecciona un nucleótido aleatorio diferente al original de forma uniforme.
    opciones = ['A', 'C', 'G', 'T']
    opciones.remove(nucleotido)      # Elimina el nucleótido original
    return random.choice(opciones) # Elige uno de los restantes

"""
Cuenta cuantas veces se repite un nucleotido en codeword recibido
"""
def freq_nucleotidos(cod_recibido):
    # Contar frecuencias de nucleótidos en recibido
    frecuencias = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '*': 0}
    #Recorre la cadena de nucleotidos, para cada nucleotido se incrementa en +1 en su correspondiente categoría 
    for nuc in cod_recibido:
        frecuencias[nuc] += 1
    return frecuencias 

"""
Decodifica usando majority voting, teniendo en cuenta las frecuencias de cada nucleotido
"""

def fill_erasures(secuencia,prior):
    secuencia_filled = []
    for nuc in secuencia:
        if nuc == '*':
            # Elegir nucleótido aleatorio según prior
            nuevo_nuc = rand_nucleotido(prior)
            secuencia_filled.append(nuevo_nuc)
        else:
            secuencia_filled.append(nuc)
    return secuencia_filled
    
def majority_voting(nucleotido, recibido,prior):
    validos = fill_erasures(recibido,prior)
    freq = freq_nucleotidos(validos)
    #print(freq)
    #print(f"maxxxxx{max(freq.values())}")     
    if max(freq.values()) < (len(validos)/2): #Si no hay candidato que se repita minimo 2 veces, aleatorio
        res = list(dict.fromkeys(validos))
        probabilidades_relativas = []
        for n in res:
            probabilidades_relativas.append(prior[n])
        #print(res)
        #print(f"priorssss{probabilidades_relativas}")
        #print(f"maximos a elegir entre{recibido}")
        return random.choices(res,probabilidades_relativas)[0] 
        
    else: #Si existe un nuc que se repite 2 ó 3 veces se elige
        return max(freq, key=freq.get)
"""
Calcula la probabilidad de transición dado un nucleotido
"""
def prob_transicion(nuc,recibido,p,n):
    #Cuento cuantos nucleotidos correctos tengo en mi secuencia
    conteo_nuc = recibido.count(nuc)

    # Si todos los nucleotidos son correctos AAA
    if conteo_nuc == n: 
        return((1-p)**conteo_nuc)
    # Si la secuencia tiene algun nucleotido correcto e incorrecto AXX AAX
    if conteo_nuc != 0 and conteo_nuc != n: 
        return((((1-p)**conteo_nuc))*((p/3)**(n-conteo_nuc)))
    # Si todos los nucleotidos estan incorrectos XXX mn    
    if conteo_nuc == 0:
        return((p/3)**n)

"""
Devuelve la log(prob_transicion) dado un nucleotido y la secuencia recibida después de pasar por el canal 
"""
def log_prob_transicion(nuc,recibido,p,n):
    epsilon = 1e-80
    # Escoge el maximo entre la probabilidad y un epsilon para evitar el cero
    log1 = math.log(max(1 - p, epsilon))
    log2 = math.log(max(p / 3, epsilon))

    conteo_nuc = recibido.count(nuc)

    return conteo_nuc * log1 + (n - conteo_nuc) * log2

"""
Dados los nucleotidos maximos, elige al azar en caso de empate 
"""
def choose_maximum(max_nuc_candidates,decoder,prior):
    probabilidades_relativas = []
    for n in max_nuc_candidates:
        probabilidades_relativas.append(prior[n])
        
    if len(max_nuc_candidates) > 1:
        if decoder == 'ML' or decoder == 'logML'or decoder == 'logMAP':
            return random.choice(max_nuc_candidates)[0]
        elif decoder in ('MAP'):
            return random.choices(max_nuc_candidates,weights = probabilidades_relativas)[0]
            
    elif len(max_nuc_candidates) == 1:
        return max_nuc_candidates[0]
"""
Dada la probabilidad de transición para cada nucleotido, se detecta cual es el maximo 
y devuelve los nucleotidos maximos
"""
def max_nucleotide(probabilidades,decoder):
    #print(f"probs_max_nucleotide {probabilidades}")
    max_nucleotide = [] 
    max_prob = max(probabilidades.values()) # Detectamos el valor maximo en probabilidad
    eps=1e-9
    if decoder in ("logML","logMAP"):
        for nuc, val in probabilidades.items():
            if abs(val - max_prob) < eps:      # ¿Está prácticamente igual al máximo?
                max_nucleotide.append(nuc)
        #print(f"max nuc logML/MAP {max_nucleotide}")
        return max_nucleotide
    else:
        for nucleotide in probabilidades:# Seleccionamos los valores que tengan ese maximo
            #print(f"nuclearizaríais:{nucleotide}, problema {probabilidades[nucleotide]}")
            if probabilidades[nucleotide]!= max_prob:
                continue
            elif probabilidades[nucleotide]== max_prob:
                max_nucleotide.append(nucleotide)
        return max_nucleotide
"""
Calcula la probabilidad de transición para cada nucleotido
"""
def all_nuc_probs(prior, recibido, decoder,p,n):
    probabilidades = {}
    for elem  in ['A', 'C', 'G', 'T']:
        if decoder == 'MAP':
            probabilidades[elem] = prior[elem] * prob_transicion(elem, recibido,p,n)

        elif decoder == 'ML':
            probabilidades[elem] = prob_transicion(elem, recibido,p,n)

        elif decoder == 'logML':
            probabilidades[elem] = log_prob_transicion(elem, recibido, p, n)
            #print(f"log_prob: {elem} { probabilidades[elem]} ")
        elif decoder == 'logMAP':
            log_prior = math.log(prior[elem])
            log_likelihood = log_prob_transicion(elem, recibido, p, n)
            probabilidades[elem] = log_prior + log_likelihood
    #print(probabilidades)
    return probabilidades
"""
Calcula las probabilidades de transición para cada nucleotido, detecta los nucleotidos con maxima probabilidad
y elige el maximo al azar si es necesario
"""
def prob_decoders(prior, recibido, decoder,p,n):

    validos = fill_erasures(recibido,prior)

    #obtener probabilidades de transición para cada nucleotido
    nucs_probs = all_nuc_probs(prior, validos, decoder,p,n)
    #Dadas las probabilidades, devuelve que nucleotidos maximizan la probabilidad
    max_nucleotides = max_nucleotide(nucs_probs,decoder)
    #Dados los candidatos se elige ya sea el unico o al azar entre los empates
    decoded_nuc = choose_maximum(max_nucleotides,decoder,prior)
    return decoded_nuc
"""
Simula la transmision sobre un canal simetrico con probabilidad p de cambio de nucleotido
"""
def channel_transmision(enviado,original,p):
    recibido = list(enviado)
    for nuc in range(len(enviado)):
        if rand_prob() > p: #Si no hay cambio de letra, continuamos al siguiente nucleotido
            continue
        else:
            recibido[nuc] = rand_nuc_diferente(original)
    return recibido
"""
Simula la transmision sobre un canal simetrico con probabilidad p de cambio de nucleotido
"""
def bio_channel(enviado,original,p_subs,p_erase,p_ins):
    #print(f"Bio channel P_subs:{p_subs} P_erase:{p_erase} P_ins:{p_ins} ")
    recibido = []
    p_max = 1 - p_erase - p_ins
    if p_subs>p_max:
        p_subs=p_max
    for nuc in range(len(enviado)):
        u = rand_prob()
        if u < p_erase: # erase
            recibido.append('*')
            #print("erase")
            
        elif u < p_erase + p_ins: # insert
            recibido.append(enviado[nuc])
            recibido.append(random.choice(['A','C','G','T']))
            #print("insert")    
        elif u < p_erase + p_ins + p_subs:#Substitute
            recibido.append(rand_nuc_diferente(original))
            #print("subs")
        else:
            recibido.append(enviado[nuc])
            #print("ok")
    return recibido
def ins_channel(enviado,original,p_subs,p_ins):
    #print(f"Bio channel P_subs:{p_subs} P_erase:{p_erase} P_ins:{p_ins} ")
    recibido = []
    for nuc in range(len(enviado)):
        u = rand_prob()
        if u < p_ins: 
            recibido.append(enviado[nuc])
            recibido.append(random.choice(['A','C','G','T']))
    
        elif u < p_ins + p_subs:#Substitute
            recibido.append(rand_nuc_diferente(original))
            #print("subs")
        else:
            recibido.append(enviado[nuc])
            #print("ok")
    return recibido

def erase_channel(enviado,original,p_subs,p_erase):
    #print(f"Bio channel P_subs:{p_subs} P_erase:{p_erase} P_ins:{p_ins} ")
    recibido = []
    for nuc in range(len(enviado)):
        u = rand_prob()
        if u < p_erase: # erase
            recibido.append('*')
            #print("erase")
  
        elif u < p_erase + p_subs:#Substitute
            recibido.append(rand_nuc_diferente(original))
            #print("subs")
        else:
            recibido.append(enviado[nuc])
            #print("ok")
    return recibido

"""
Retorna 1 en caso de error
"""
def contar_error(original, decodificado):
    return 1 if original != decodificado else 0

######################################################################################################################################################
puntos = 20
#p_vals = p_curve = np.logspace(-3, 0, puntos)  # 50 puntos entre 10^-3 y 1
#p_vals = p_curve = np.linspace(0.001, 1, puntos)  # Por ejemplo
p_vals = p_curve = np.linspace(0.001, 1, puntos)  # Por ejemplo

#p_curve = np.linspace(0.6, 1, puntos)  # Por ejemplo
#p_vals =[0.75]
sim_mv_erase,sim_mv_ins,sim_mv,sim_mv_ins, sim_ml, sim_map, sim_logml,sim_logmap = [], [], [], [], [], [], [], []
sim_mv_bio, sim_ml_bio, sim_map_bio, sim_log_ml_bio, sim_log_map_bio = [], [], [], [], []
contador = 0
p_erase = 0
p_ins = 0
n = 1 #Longitud de codigo de repetición
prior = {
    'A': 0.1,
    'C': 0.4,
    'G': 0.4,
    'T': 0.1
}

for p in p_vals: # Barrido de p
    print(datetime.now().strftime("%H:%M:%S"))
    #N = N_requerido(p, k=10)     # usa k=100 o 1000 según precisión deseada
    N = 200_000
    print(f"p = {p:.1e}   →   Pe≈{Pe_teorico(p):.1e}   →   N = {N:.1e}")
    
    #print(p)
    print(f"progreso:{contador}/{puntos}")
    contador+=1
    error_mv_erase = error_mv_ins =  error_mv = error_mv_ins = error_ml = error_map = error_logml = error_logmap = 0
    
    error_mv_bio = error_ml_bio = error_map_bio= error_log_ml_bio= error_log_map_bio = 0
    
    for _ in range (N): # Para N iteraciones 
        nucleotido = rand_nucleotido(prior) #Creamos el nucleotido de manera aleatoria
        deco_nucleotido = ' ' # variable para almacenar el nucleotido decodificado
    
#________________Transmision
        sent_codeword = nucleotido * n # codigo de repetición
        """
        received_codeword = channel_transmision(sent_codeword,nucleotido,p) # Efectos del canal
        received_erase = erase_channel(sent_codeword,nucleotido,p,p_erase)
        received_ins = ins_channel(sent_codeword,nucleotido,p,p_ins)
        """
        #received_bio = bio_channel(sent_codeword,nucleotido,p,p_erase,p_ins)
        received_bio = bio_channel(sent_codeword,nucleotido,p,p_erase,p_ins)
        
        #print(p)
       # print(f" enviado, recibido {sent_codeword},{received_bio} ")

        
#_________________Decodificación 
        """
        deco_mv_erase = majority_voting(nucleotido, received_erase,prior)# ok 
        deco_mv_ins = majority_voting(nucleotido, received_ins,prior)
        
        deco_mv =  majority_voting(nucleotido, received_codeword,prior)
        deco_ml = prob_decoders(prior, received_codeword, 'ML',p,n)
        deco_map = prob_decoders(prior, received_codeword, 'MAP',p,n)
        deco_logml = prob_decoders(prior, received_codeword, 'logML', p, n)
        deco_logmap = prob_decoders(prior, received_codeword, 'logMAP', p, n)
        """

        
        deco_mv_bio = majority_voting(nucleotido, received_bio,prior)
        """
        deco_ml_bio = prob_decoders(prior, received_bio, 'ML',p,n)
        deco_map_bio = prob_decoders(prior, received_bio, 'MAP',p,n)
        """
        deco_log_ml_bio = prob_decoders(prior, received_bio, 'logML',p,n)
        deco_log_map_bio = prob_decoders(prior, received_bio, 'logMAP',p,n)

#_________________Comprobacion de Decodificación
        error_mv_bio +=contar_error(nucleotido,deco_mv_bio) 
        """
        error_ml_bio +=contar_error(nucleotido,deco_ml_bio) 
        error_map_bio +=contar_error(nucleotido,deco_map_bio) 
        """
        error_log_ml_bio +=contar_error(nucleotido,deco_log_ml_bio) 
        error_log_map_bio +=contar_error(nucleotido,deco_log_map_bio) 
        
        
        """
        error_mv_erase +=contar_error(nucleotido,deco_mv_erase)
        error_mv_ins +=contar_error(nucleotido,deco_mv_ins)

        
        error_mv += contar_error(nucleotido,deco_mv)
        error_ml += contar_error(nucleotido,deco_ml)
        error_map += contar_error(nucleotido,deco_map)
        error_logml += contar_error(nucleotido, deco_logml)
        error_logmap += contar_error(nucleotido, deco_logmap)

        """
    # Guardamos promedios
    sim_mv_bio.append(error_mv_bio/N)
    sim_ml_bio.append(error_ml_bio/N)
    sim_map_bio.append(error_map_bio/N)
    sim_log_ml_bio.append(error_log_ml_bio/N)
    sim_log_map_bio.append(error_log_map_bio/N)
"""


    sim_mv_erase.append(error_mv_erase / N)
    sim_mv_ins.append(error_mv_ins / N)
    
    sim_mv.append(error_mv / N)
    sim_ml.append(error_ml / N)
    sim_map.append(error_map / N)
    sim_logml.append(error_logml / N)
    sim_logmap.append(error_logmap / N)

"""
#------------------------------------------------------------------------------------------------------------


Error_Teorico_MV = (7 * p**2 - 4 * p**3) / 3   



# Crear la figura
plt.figure(figsize=(10, 6))
ax = plt.gca()
ax.set_facecolor('white')
# Curvas teóricas

theo_mv = (7 * p_curve**2 - 4 * p_curve**3) / 3
theo_ml = np.where(p_curve < 0.75, theo_mv, 1 - 16 * (p_curve / 3)**3)

plt.plot(p_curve, theo_mv, label='Theoretical Substitution Model (MV Decoder)', color='#485460', linewidth=1.5, linestyle='--')
plt.plot(p_curve, theo_ml, label='Theoretical Substitution Model (ML Decoder)', color='#1e272e', linewidth=1.5)

# Puntos simulados
plt.scatter(p_vals, sim_mv_bio,
            label='MV ',
            color='#05c46b',
            marker='.', s=40)
plt.plot(p_vals, sim_mv_bio,
         color='#05c46b',
         linestyle='-',
         linewidth=1)
"""
plt.scatter(p_vals, sim_ml_bio,
            label='ML ',
            color='#ffa801',
            marker='x', s=40)
plt.plot(p_vals, sim_ml_bio,
         color='#ffa801',
         linestyle='-',
         linewidth=1)

plt.scatter(p_vals, sim_map_bio,
            label='MAP ',
            color='#f53b57',
            marker='_', s=40)
plt.plot(p_vals, sim_map_bio,
         color='#f53b57',
         linestyle='-',
         linewidth=1)

"""
plt.scatter(p_vals, sim_log_ml_bio,
            label='Log ML ',
            color='#D980FA',
            marker='_', s=40)
plt.plot(p_vals, sim_log_ml_bio,
         color='#D980FA',
         linestyle='-',
         linewidth=1)
plt.scatter(p_vals, sim_log_map_bio,
            label='Log MAP',
            color='#0fbcf9',
            marker='x', s=40)
plt.plot(p_vals, sim_log_map_bio,
         color='#0fbcf9',
         linestyle='-',
         linewidth=1)

 


# Línea vertical en p = 0.75
plt.axvline(x=0.75, color='gray', linestyle=':', linewidth=2, label='p = 0.75')

# Escalas logarítmicas
#plt.xscale('log')
#plt.yscale('log')

# Etiquetas y leyenda
plt.xlabel('p')
plt.ylabel('Error Probability')
plt.title(f'BioChannel \n erasure = {p_erase*100:.0f}%, insert = {p_ins*100:.0f}%, Repetition Code (Length = {n})',fontsize=12, color='#1e272e')
bbox=dict(facecolor='white', alpha=0.8)
plt.grid(True, which='both', linestyle=':')
plt.legend()
plt.tight_layout()

plt.gca().set_aspect('auto')  # relación aspecto automática
plt.savefig('bio.pdf', format='pdf', bbox_inches='tight' )
plt.show()
  