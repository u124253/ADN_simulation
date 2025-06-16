# Decodificación en ADN
Este proyecto genera simulaciónes sobre modelos de almacenamiento biologico en cadenas de ADN utilizando codigos de repetición y distintos decodificadores.
## Modelos:
-  Canal de sustituciónes
-  Canal de borrados
-  Canal de inserciónes
-  Canal biológico (combinacion de los 3 anteriores)

## Requisitos
	•	Python 3.x
	•	matplotlib

Instalación rápida de dependencias:
pip install matplotlib

## Uso
Ejecuta el script principal: python BioChannel.py

Los archivos PDF generados se guardarán automáticamente con los parametros seleccionados:
CanalBiolofico_n100_ins10_borrados10.pdf

## Parámetros modificables:
  -  n: longitud del código de repetición 
  -  p_ins: probabilidad de inserción (valor entre 0 y 1)
  -  p_erase: probabilidad de borrado (valor entre 0 y 1)
  -  N , K: número de interaciones deseadas y su presición 

## Decodificadores:
 -  Majority Voting
 -  Maximum Likelihood
 -  Maximum a Posteriori
### Para la estimación de probabilidades muy pequeñas se puede seleccionar:
       - Maximum Likelihood logaritmico
       - Maximum a Posteriori Logaritmico

## Salida
Los gráficos se guardan en formato PDF en el mismo directorio que el script.

### Licencia
Este proyecto está bajo la licencia MIT.
