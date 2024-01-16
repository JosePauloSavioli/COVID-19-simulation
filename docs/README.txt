Chamada do programa 'simulacao.py'

---> Instalar as bibliotecas necessárias (ver arquivo "requirements.txt")
---> Ir para o diretório da simulação onde está o arquivo 'simulacao.py'
---> Checar se há uma pasta chamada 'Simulações' no diretório raiz

---> Sintaxe:

'python' ou 'python3'
+
'simulacao.py'
+
'-'código do município (no caso de análises em loop pelo csv este campo pode ser deixado vazio ['-']) [ex: -28]
+
'-'nome do arquivo com os casos confirmados (sem '.csv') [ex (conf.csv): -conf]
+
'-'nome do arquivo com as mortes acumuladas (sem '.csv') [ex (mort.csv): -mort]
+
'-'método de fitting [-0: basinhopp, -1: differential evolution [default], -2: powell, -3: cobyla] [ex: -1]
+
'-'booleana para utilização de tau fantasma [-0: não utilizar, -1: utilizar]
'-'quantidade de taus fantasmas [ex: -0-0 (sem tau fantasma), -1-2 (utilização de 2 taus fantasmas no fitting)]
'-'booleana para utilização de obt fantasma [-0: não utilizar, -1: utilizar]
'-'quantidade de obts fantasmas [ex: -0-0 (sem tau fantasma), -1-2 (utilização de 2 taus fantasmas no fitting)]
+
'-'tipo de simulação [-n: simulação de uma localidade, -s: simulação de todas as localidades de um .csv, -b: bootstrap de uma localidade, -sl: simulação de uma localidade com sensibilidade] [ex: -n]
+
'-'período de simulação em dias [ex: -200]
+
'-'número de dias para validação [ex: -5]
+
'-'tipo de simulação [-mod: simulação com hospitalizados, -std: simulação SEIR com assintomáticos e mortes]
+
'-'rodar gráficos e testes adicionais [-0: não, -1: sim]

---> Exemplo de chamada (casos e óbitos nos estados do Brasil)
python simulacao.py -28 -nc85 -oa85 -1 -1-2-0-0 -b -200 -5 -str -0

