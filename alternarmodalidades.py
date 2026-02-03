# version simple, testear primero esta.


# mejor version: usando un generador

import os

def mi_funcion_pesada(sublista, indice):
    nombre_carpeta = f"resultado_sin_elemento_{indice}"
    os.makedirs(nombre_carpeta, exist_ok=True)

    # Supongamos que guardamos un archivo dentro
    with open(f"{nombre_carpeta}/data.txt", "w") as f:
        f.write(f"Procesado con {len(sublista)} elementos")

    # No retornamos nada pesado, solo confirmación
    return f"Carpeta {nombre_carpeta} creada."

elementos = ["A", "B", "C", "D"]

# El generador asegura que no satures la RAM
for i, sublista in enumerate(elementos[:j] + elementos[j+1:] for j in range(len(elementos))):
    status = mi_funcion_pesada(sublista, i)
    print(status) # Solo para ver el progreso


# correr version completa
modalidadesorg = ["Metabolomics", "Lipidomics", "Microbiota", "SV", "VC", "AEB", "AS"]
main(modalidadesorg, "CompleteRun")

# correr quitando una modalidad
modalidades = ["Metabolomics", "Lipidomics", "Microbiota", "SV", "VC", "AEB", "AS"]

primervalor = modalidades[0]
fr = modalidades.pop()
modalidades.insert(0, fr)

while primervalor != modalidades[0]:
    # quita el último elemento
    fr = modalidades.pop()
    # da un nombre a la run
    withoutmod = list(set(modalidadesorg) - set(modalidades))
    withoutmod = ", ".join(map(str, withoutmod))
    runname = f"runwithout{withoutmod}"
    # ejecuta el script tocho con las modalidades del txt
    main(modalidades, runname)
    # añade el elemento descartado a la primera posición
    modalidades.insert(0, fr)


# lo suyo sería que el nombre de la run fuera "runwithoutx" siendo x la modalidad que no está en la lista
# para ello al correr main se le tiene que poder dar un nombre como variable que saldra de aquí
withoutmod = list(set(modalidadesorg) - set(modalidades))
withoutmod = ", ".join(map(str, withoutmod))
runname = f"runwithout{withoutmod}"

# ademas las guiding variables van a cambiar, porque tienes varios subgrupos
#We have a lot of metadata that can be our guiding variables. The idea is to separate our metadata in clinic variables that Sergi pointed out as clinically related to EOCRC. That defines 6 different groups. Therefore, we'll do 6 different analysis per try. 
#The groups are:
# Lifestyle factors}: tabaquism and alcoholism.
# Tumor localization}: según si es ascendent, transversal or recto.
# Tumor characteristics}: grade of diferentiation and mucinous content.
# Cancer stage}: stages I, II, III or IV.
# Metastasis}: presence or absence.
# Medical record}: presence of previous polips and CRC in the family.


# version variando factores
modalidadesorg = ["Metabolomics", "Lipidomics", "Microbiota", "SV", "VC", "AEB", "AS"]
modalidades = ["Metabolomics", "Lipidomics", "Microbiota", "SV", "VC", "AEB", "AS"]
factors = [10, 20, 30]

primervalor = modalidades[0]
fr = modalidades.pop()
modalidades.insert(0, fr)

while primervalor != modalidades[0]:
    # quita el último elemento
    fr = modalidades.pop()
    # escribe modalidades en un txt
    with open("salida.txt", "w") as archivo:
        archivo.write(modalidades)
    # ejecuta el script tocho con las modalidades del txt
    for i in factors:
        main(modalidades, i)
    # añade el elemento descartado a la primera posición
    modalidades.insert(0, fr)



# version recursiva a corregir
modalidades = ["Metabolomics", "Lipidomics", "Microbiota", "SV", "VC", "AEB", "AS"]
def recmodalidades():
    primervalor = modalidades[0]
    fr = modalidades.pop()
    modalidades.insert(0, fr)

    while primervalor != modalidades[0]:
        fr = modalidades.pop()
        # escribe modalidades en un txt
        print(modalidades)
        recmodalidades()
        # ejecuta el script tocho con las modalidades del txt
        modalidades.insert(0, fr)
recmodalidades()
