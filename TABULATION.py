import pandas as pd #Se tuvo que descargar el módulo antes con pip install pandas
import os #No se descarga, viene por defecto. Permite navegar/crear/listar archivos por carpetas

#EJECUTAR EL SCRIPT CON python AL PRINCIPIO PARA QUE RECONOZCA EL DEL ENTORNO QUE ES EL QUE TIENE PANDAS INSTALADO


#lista_parental= ["ATCC49226", "M6_S40"]

ruta = os.path.expanduser('~/Proyectos/Prueba_programas_actualizados2/curated_files_ATCC49226') #Esto sirve para que reconozca ~
nombres_curated = [] #Se define la lista vacia para que lo reconozca el bucle


for archivo in os.listdir(ruta): #lista todos los archivos que hay en la ruta determinada
    if archivo.endswith("_curated") and not archivo.endswith("_pre_curated"): #Si el archivo termina con _curated y no con _pre_curated, continua
        nombre = archivo.split("_curated")[0] #split recorta el nombre por _curated, y [0] escoge lo que va antes del corte
        nombres_curated.append(nombre) #Almacenamos cada nombre en una lista
        ruta_completa = os.path.join(ruta, archivo) #os.path.join une varias partes de una ruta de archivo en una sola usando el separador /
        print(f"Leyendo {ruta_completa}")
        df = pd.read_csv(ruta_completa, sep="\t") #Convierte en df cada archivo
        #print(df.head(2))

        df["CHROM"] = df["CHROM"].str.replace("ATCC_19424_", "", regex=False)
        df["CHROM"] = df["CHROM"].apply(lambda x: f"{nombre}_{x}") #.apply aplica una función a cada fila y lambda permite hacer una función anónima a la fila (x). Así, con el f-string se "cambia" la fila por lo que hay guardado en nombre y x que es lo que habia en la fila originalmente
        print(df.tail(5))
        
print(nombres_curated)