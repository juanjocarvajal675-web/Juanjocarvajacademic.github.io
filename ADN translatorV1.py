# Transcripción, replicación, traducción y expresión del ADN
# Autor: Juan José Carvajal Marín
# Descripción: Este proyecto es una herramienta de consola escrita en Python (versión 3.10 o superior) diseñada para simular los procesos fundamentales del dogma central
# de la biología molecular a partir de secuencias genéticas en dirección 5' a 3'. Para utilizarlo, ejecuta el archivo desde tu terminal y selecciona una de las cinco operaciones
# del menú principal ingresando su número correspondiente: transcripción de ADN a ARN (opción 1), replicación de ADN a su hebra complementaria (opción 2), transcripción inversa de
# ARN a ADN (opción 3), traducción directa de ARN a proteínas (opción 4) o traducción completa desde ADN hasta proteínas (opción 5). Si eliges las opciones de traducción (4 o 5),
# el programa te pedirá elegir entre una traducción completa —que convierte toda la secuencia en bloques de tres nucleótidos sin filtros— o una traducción inteligente,
# la cual detecta automáticamente cada marco de lectura abierto iniciando exclusivamente en el codón de arranque AUG (Metionina) y deteniéndose en el primer codón de parada
# (UAA, UAG o UGA) para entregarte las proteínas resultantes tanto en formato de tres letras como en su nomenclatura abreviada de una letra. Asegúrate de ingresar únicamente las
# bases nitrogenadas válidas (A, T, C, G para ADN; A, U, C, G para ARN) sin espacios ni caracteres especiales, ya que de lo contrario el sistema detendrá la ejecución indicando que
# hay nucleótidos inválidos.
def transcripcion_ADN(cadena = None):
    cadena = str(input("Ingrese la cadena de ADN a trascribir: "))
    cadena = cadena.upper()
    transcrito = "0"
    try:
        for nucleotido in cadena:
            transcrito += Bp_ADN_ARN[nucleotido]
        return (transcrito[1:])
    except: print("Nucleotidos inválidos")
    
def transcripcion_inv_ARN(cadena = None):
    cadena = str(input("Ingrese la cadena de ARN a trascribir: "))
    cadena = cadena.upper()
    transcrito = "0"
    try:
        for nucleotido in cadena:
            transcrito += Bp_ARN_ADN[nucleotido]
        return (transcrito[1:])
    except: print("Nucleotidos inválidos")
def replicacion_ADN(cadena = None):
    cadena = str(input("Ingrese la cadena de ADN a trascribir: "))
    cadena = cadena.upper()
    replicado = "0"
    try:
        for nucleotido in cadena:
            replicado += Bp_ADN_ADN[nucleotido]
        return (replicado[1:])
    except: print("Nucleotidos inválidos")
    
def traduccion_ARN(cadena = None):
    tipo = int(input("1. traducción completa de los codones\n2. traducción inteligente\n>> "))
    if cadena == None:
        cadena = str(input("Ingrese la cadena de ARN a trascribir en direción 5´ a 3´: "))
    cadena = cadena.upper()
    
    codones = [cadena[i:i+3] for i in range(0, len(cadena) - len(cadena) % 3, 3)]
    print(codones)
    
    match tipo:
        case 1:
            traducido = ""
            for codon in codones:
                traducido += Bp_ARN_codon[codon] + " "
            resultado = traducido.strip()
            print(resultado)
            print("Resumida:", resumir_proteina(resultado))
            return None
        case 2:
            proteinas = []
            traducido = ""
            activo = False
            for codon in codones:
                if Bp_ARN_codon[codon] == "Met" and not activo:
                    activo = True
                if activo:
                    if Bp_ARN_codon[codon] == "STOP":
                        activo = False
                        proteinas.append(traducido.strip())
                        traducido = ""
                        continue
                    traducido += Bp_ARN_codon[codon] + " "
            return proteinas

def resumir_proteina(aminoacidos):
    lista_aa = aminoacidos.split()
    res_proteina = ""
    for aa in lista_aa:
        res_proteina += Codon_resumido[aa]
    return res_proteina
    
    
    
    
# Define los pares complementarios (todo está en mayúsculas)
Bp_ADN_ADN = {"T":"A","A":"T","C":"G","G":"C"}
Bp_ADN_ARN = {"T":"A","A":"U","C":"G","G":"C"}
Bp_ARN_ADN = {"U":"A","A":"T","C":"G","G":"C"}
Bp_ARN_codon = codones = {
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
    "AAU": "Asn", "AAC": "Asn",
    "GAU": "Asp", "GAC": "Asp",
    "UGU": "Cys", "UGC": "Cys",
    "CAA": "Gln", "CAG": "Gln",
    "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
    "CAU": "His", "CAC": "His",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile",
    "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AAA": "Lys", "AAG": "Lys",
    "AUG": "Met",
    "UUU": "Phe", "UUC": "Phe",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "AGU": "Ser", "AGC": "Ser",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "UGG": "Trp",
    "UAU": "Tyr", "UAC": "Tyr",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UAG": "STOP", "UGA": "STOP", "UAA": "STOP",
}
Codon_resumido = {
    "Ala": "A",   # Alanina
    "Arg": "R",   # Arginina
    "Asn": "N",   # Asparagina
    "Asp": "D",   # Ácido aspártico
    "Cys": "C",   # Cisteína
    "Gln": "Q",   # Glutamina
    "Glu": "E",   # Ácido glutámico
    "Gly": "G",   # Glicina
    "His": "H",   # Histidina
    "Ile": "I",   # Isoleucina
    "Leu": "L",   # Leucina
    "Lys": "K",   # Lisina
    "Met": "M",   # Metionina (inicio)
    "Phe": "F",   # Fenilalanina
    "Pro": "P",   # Prolina
    "Ser": "S",   # Serina
    "Thr": "T",   # Treonina
    "Trp": "W",   # Triptófano
    "Tyr": "Y",   # Tirosina
    "Val": "V",   # Valina
    "STOP": "*",  # Codón de parada
}


# selector de operación
try:
    opcion = int(input("1. trancripción ADN a ARN\n2. replicación de ADN a ADN\n3. transcrición inversa de ARN a ADN\n4. tradución de ARN a proteínas (aminoácidos)\n5. tradución de ADN a proteínas (aminoácidos)\n>> "))
    match opcion:
        case 1: print(transcripcion_ADN())
        case 2: print(replicacion_ADN())
        case 3: print(transcripcion_inv_ARN())
        case 4:
            proteinas = traduccion_ARN()
            if proteinas:
                for i in range(len(proteinas)):
                    print(f"Proteína {i+1} ({resumir_proteina(proteinas[i])}): {proteinas[i]}")
        case 5:
            proteinas = traduccion_ARN(transcripcion_ADN())
            if proteinas:
                for i in range(len(proteinas)):
                    print(f"Proteína {i+1} ({resumir_proteina(proteinas[i])}): {proteinas[i]}")

        case _: print("Opción incorrecta")
except ValueError:
    print("Opción incorrecta")
except KeyError:
    print("Nucleótidos inválidos en la cadena")
        
