## Build kaiju database ##

mamba activate kaiju

kaiju-mkbwt -a ACDEFGHIKLMNPQRSTVWY -o proteins $proteins_kaiju
kaiju-mkfmi proteins
