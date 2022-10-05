# Modelo de Ising 2D - Serial
[![NPM](https://img.shields.io/npm/l/react)](https://github.com/DanSantiago/Modelo-de-Ising-2D/blob/main/LICENCE) 

# Sobre o projeto

O projeto foi desenvolvido durante a disciplina de computação de alto desempenho I para desenvolver os trabalhos durante o período. O algoritmo consiste em simular
o modelo de Ising 2D em serial.

## Composição 
 
O projeto é composto por três arquivos: "funcoes.c", "funcoes.h" e "ising.c".

- funcoes.h: header
- funcoes.c: biblioteca com as funções utilizadas
- ising.c: programa principal

As constantes e a inclusão de bibliotecas utilizadas no modelo são definidas/incluídas no header (funcoes.h) para os outros arquivos terem acesso a elas.

# Tecnologias utilizadas
- C

## Bibliotecas utilizadas
- stdio
- stdlib
- math

OBS: O algoritmo foi gerado no ambiente Linux e podem ocorrer erros de compatibilidade se executado no ambiente Windows.

# Como executar o projeto

Após baixar a pasta do projeto, basta seguir os seguintes passos:

Inicialmente é necessário configurar no arquivo makefile o compilador ao qual deseja utilizar e as flags de interesse.

No compilador gcc há uma flag básica: -lm.

No compilador ifort não há necessidade de nenhuma flag básica.

As flags devem ser definidas no campo CFLAGS no arquivo makefile.

Para usar o compilador da intel defina icc no campo CC do arquivo makefile.

Para usar o compilador da gnu defina gcc no campo CC do arquivo makefile.

Depois da configuração ter sido realizada, basta utilizar o comando "make" e diversos arquivos serão gerados, inclusive o executável ising.x.

Após isso, basta executar o arquivo criado e o modelo irá rodar, salvando os valores no arquivo dados.dat na seguinte ordem: Temperatura | Magnetização | Energia | Calor específico | Susceptibiliade Magnética.

Após a execução, basta utilizar o comando "make clean" para apagar os arquivos gerados para a execução do programa.

Algumas variáveis de execução podem ser modificadas caso seja de interesse, uma breve descrição delas:

- L: tamanho da malha
- N: número de spins
- Nskip: passos para termalização
- T: temperatura máxima
- Tmin: temperatura mínima
- dT: incremento da temperatura
- Msteps: número de passos de Monte Carlo
- J: valor da constante de troca
- k: valor da constante de Boltzmann
- B: valor do campo magnética externo

# Autor

Daniel Santiago da Silva

https://www.linkedin.com/in/danielssilva13/
