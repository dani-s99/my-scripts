import pandas as pd #Se tuvo que descargar el módulo antes con pip install pandas
import os #No se descarga, viene por defecto. Permite navegar/crear/listar archivos por carpetas

#EJECUTAR EL SCRIPT CON python AL PRINCIPIO PARA QUE RECONOZCA EL DEL ENTORNO QUE ES EL QUE TIENE PANDAS INSTALADO

ruta_principal = os.path.expanduser('~/Documentos/Examples/Subflava')#Esto sirve para que reconozca ~
nombres_curated = [] #Se definen las listas vacias para que lo reconozca el bucle
lista_df = []
parental_only = []

for carpeta in os.listdir(ruta_principal): #lista todas las carpetas que hay en la ruta determinada
    ruta_carpeta = os.path.join(ruta_principal, carpeta) # os.path.join une varias partes de una ruta de archivo en una sola usando el separador /
    if carpeta.startswith("curated_files_"): #Si la carpeta empieza por curated_file_, continua
        print (f"Se ha podido entrar a la carpeta {carpeta}")
        for archivo in os.listdir(ruta_carpeta):
            if archivo.endswith("_curated") and not archivo.endswith("_pre_curated"): #Si el archivo termina con _curated y no con _pre_curated, continua
                nombre = archivo.split("_curated")[0] #split recorta el nombre por _curated, y [0] escoge lo que va antes del corte
                nombres_curated.append(nombre) #Almacenamos cada nombre en una lista
                ruta_completa = os.path.join(ruta_carpeta, archivo)
                print(f"Leyendo {ruta_completa}")
                df = pd.read_csv(ruta_completa, sep="\t") #Convierte en df cada archivo
                #print(df.head(2))

                df["CHROM"] = df["CHROM"].apply(lambda x: f"{nombre}") #.apply aplica una función a cada fila y lambda permite hacer una función anónima a la fila (x). Así, con el f-string se "cambia" la fila por lo que hay guardado en nombre y x que es lo que habia en la fila originalmente ({nombre}_{x}). He quitado {nombre}_{x} y solo he dejado {nombre} para que no diferencie entre contig_1 y contig_2 ya que los genes que están en contig_2 solo aparecen ahí.
                print(df.tail(5))

                lista_df.append(df) #Añadir dataframe a la lista, esto es una lista, no es el archivo final con todas las filas juntas. Para esto tenemos que salir del bucle 

df_total = pd.concat(lista_df, ignore_index=True)  #Concatenar todos los df uno debajo del otro. Con ignore_index=True se ignoran los índices de cada fila que se generan por cada dataframe

# Mostrar una parte para verificar
# print("DataFrame combinado:")
# print(df_total)

# Mostrar los nombres procesados
# print("Archivos procesados:", nombres_curated)

metadata = pd.read_csv("~/Documentos/METADATA", sep="\t") #Pasamos el archivo a una variable llamada metadata
# print(metadata.head())
# print(parental_only)


df_total["Samples"] = df_total["CHROM"] #Cambio el nombre de la columna CHROM por Samples.
df_total["Samples"] = df_total["Samples"].str.split("_").str[0] #Esto hace que en la columna Samples se separe la parte del nombre de la muestra (MX) de SX.
df_total["mutation_without_c_or_n"] = df_total["ANN[*].HGVS_C"].str.replace("^(c\.|n\.)", "", regex=True)
df_total["mutation_without_c_or_n"] = (
        df_total["mutation_without_c_or_n"]
    .fillna("") # reemplaza NaN por cadena vacía
    .astype(str))# convierte todo a string ya que si no da probelmas de float

#Con .merge se combinan datos de acuerdo a las opciones/claves que pongamos. En este caso, se miran solo las columnas Muestras y Info de metadata y si hay coincidencia entre Muestras (metadata) y Samples (df_total) pasa lo que hay de Info en metadata a df_total.
df_total = df_total.merge(metadata[["Muestras", "Info"]],
                          how="left", #En que está en left es el df principal
                          left_on="Samples", #En left la columna clave es Samples
                          right_on="Muestras") #En rigth (metadata), la columna clave es Muestras
                       
print("Comprobando que se añadió la columna Info al data frame")
print(df_total.head(2))

#Añadimos los resultados finales a una tabla conjunta donde las filas son cada muestra y las columanas los ID de los locus. En cada casilla tiene que haber las mutaciones encontradas para ese locus.
tabla_final = df_total.pivot_table( 
    index=["Samples", "Info"], # las filas, como hemos puesto dos columnas el index se convierte en una tupla
    columns="ANN[*].GENE",     # las columnas
    values="mutation_without_c_or_n",  # Las mutaciones encontradas para esa muestra y locus.
    aggfunc=lambda x: ', '.join(x)  # cómo combinar si hay varias mutaciones iguales
)

#Queremos que el orden de las filas sea numérico y no alfabético. Para ello, usamos el comando sorted.
tabla_final = tabla_final.loc[
    sorted(tabla_final.index, key=lambda x: int(x[0][1:])) #x es una tupla, como ('M1', 'X'); x[0] accede al primer elemento "M1"; x[0][1:] quita la letra "M"; int(...) convierte el número a entero para ordenar correctamente.
]

print(tabla_final.head(5))

tabla_final.to_excel("tabla_final_subflava.xlsx") #CAMBIAR EL NOMBRE SEGÚN LA MUESTRA



