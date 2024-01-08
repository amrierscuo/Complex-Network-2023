% tnode_centralities.m
% Questo script analizza le centralità in una rete complessa.
% Gli indici di centralità considerati sono:

% 1. Centralità di Grado (Degree Centrality): 
%    La centralità di grado è una misura fondamentale in analisi di rete che 
%    quantifica l'importanza di un nodo basandosi sul numero di collegamenti (o archi)
%    che ha con altri nodi nella rete. Questa misura è particolarmente significativa
%    in reti non orientate, dove ciascun collegamento ha la stessa importanza.

%    In reti orientate, la centralità di grado può essere ulteriormente suddivisa
%    in 'centralità di grado in entrata' (in-degree) e 'centralità di grado in uscita' 
%    (out-degree). Il grado in entrata conta il numero di archi diretti verso il nodo, 
%    mentre il grado in uscita conta il numero di archi che partono dal nodo.

%    Un nodo con un alto grado in entrata è spesso visto come un nodo influente o 
%    popolare all'interno della rete, mentre un nodo con un alto grado in uscita può 
%    essere considerato come un grande diffusore o emittente di informazioni.
%
%    Calcolo:
%    La centralità di grado per un nodo in una rete è calcolata come il conteggio 
%    del numero totale di archi (o connessioni) che il nodo ha con altri nodi. 
%    Questo conteggio include sia i collegamenti in entrata che quelli in uscita in 
%    reti orientate, o semplicemente tutti i collegamenti in reti non orientate.

%    In MATLAB, la centralità di grado è calcolata utilizzando la funzione 'centrality'
%    applicata a un oggetto grafo (tipo 'graph' o 'digraph'). La funzione 'centrality'
%    con l'opzione 'degree' esegue il conteggio degli archi per ciascun nodo nel grafo.

%    Algoritmo MATLAB:
%    - Per un grafo non diretto ('graph'), MATLAB conta il numero di archi collegati 
%      a ciascun nodo. Questo è equivalente al conteggio del numero di nodi adiacenti.
%    - Per un grafo diretto ('digraph'), MATLAB distingue tra 'in-degree' e 'out-degree'.
%      - 'In-degree' conta il numero di archi entranti in un nodo.
%      - 'Out-degree' conta il numero di archi uscenti da un nodo.
%    - La funzione opera iterando su tutti i nodi del grafo e sommando il numero di 
%      archi per ciascun nodo. Questo processo è ottimizzato per essere efficiente 
%      anche su grandi reti.

%    Formule:
%    - Per un grafo non diretto: Degree(nodo) = Numero di archi connessi al nodo
%    - Per un grafo diretto: 
%      - In-degree(nodo) = Numero di archi che entrano nel nodo
%      - Out-degree(nodo) = Numero di archi che escono dal nodo

%    Esempio:
%    Se un nodo è connesso a 5 altri nodi in una rete, la sua centralità di grado è 5.
%    Questo significa che il nodo ha cinque connessioni dirette con altri nodi nella rete.

%    Limitazioni:
%    - La centralità di grado non tiene conto della qualità o dell'importanza dei nodi 
%      a cui un dato nodo è connesso(nessun peso).
%    - In reti grandi e complesse, un alto grado potrebbe non necessariamente
%      implicare un'alta importanza a livello globale della rete.

% 2. Centralità di Vicinanza (Closeness Centrality):
%    La centralità di Closeness è una misura che indica quanto un nodo è 'vicino' 
%    a tutti gli altri nodi in una rete. In termini più tecnici, è basata sull'inverso 
%    della somma delle distanze più corte da un nodo a tutti gli altri nodi nella rete. 
%    Un valore elevato di centralità di Closeness implica che un nodo è rapidamente 
%    raggiungibile da altri nodi, giocando un ruolo cruciale nella diffusione delle 
%    informazioni o risorse nella rete.

%    Calcolo:
%    Per calcolare la centralità di Closeness di un nodo, si sommano le distanze 
%    più corte dal nodo a tutti gli altri nodi nella rete e si prende l'inverso 
%    di questa somma. In formule, se d(i, j) è la distanza più breve dal nodo i 
%    al nodo j, allora la centralità di vicinanza C(i) è data da:
%    C(i) = 1 / (Σ d(i, j) per tutti j ≠ i)

%    Algoritmo MATLAB:
%    - MATLAB calcola la centralità di Closeness utilizzando la funzione 'centrality'
%      con l'opzione 'closeness' su un oggetto grafo.
%    - Internamente, MATLAB utilizza algoritmi di ricerca del cammino minimo, 
%      come l'algoritmo di Floyd-Warshall o Dijkstra, per calcolare le distanze più 
%      brevi tra i nodi.
%    - La somma di queste distanze è poi utilizzata per calcolare l'inverso della 
%      somma totale, determinando la centralità di Closeness per ciascun nodo.

%    Limitazioni:
%    - In reti sparse o in reti con molti cluster separati, la centralità di Closeness 
%      può non essere significativa, poiché alcune parti della rete potrebbero essere 
%      irraggiungibili da certi nodi.
%    - La centralità di vicinanza assume che tutte le connessioni siano ugualmente 
%      importanti, il che potrebbe non essere vero in tutte le applicazioni.

% 3. Centralità di Intermediazione (Betweenness Centrality):
%    La centralità di Betweenness è una misura cruciale in analisi di rete 
%    che quantifica il ruolo di un nodo come intermediario nei percorsi più brevi 
%    tra coppie di nodi nella rete. Un alto valore di centralità di Betweenness 
%    indica che un nodo agisce frequentemente come punto di passaggio (o 'ponte') 
%    in molti percorsi più brevi, rendendolo strategico per il controllo del flusso 
%    di informazioni o risorse nella rete.

%    Calcolo:
%    Per calcolare la centralità di Betweenness di un nodo, si considera il numero 
%    di volte che il nodo si trova sui percorsi più brevi tra tutte le coppie di nodi 
%    nella rete. La formula per calcolarla è data da:
%    C_B(i) = Σ (σ(j, k | i) / σ(j, k)) per tutte le coppie j, k dove j ≠ k ≠ i
%    dove σ(j, k) è il numero totale di percorsi più brevi dal nodo j al nodo k 
%    e σ(j, k | i) è il numero di tali percorsi che passano attraverso il nodo i.

%    Algoritmo MATLAB:
%    - MATLAB calcola la centralità di Betweennes usando la funzione 'centrality'
%      con l'opzione 'betweenness' su un oggetto grafo.
%    - Internamente, MATLAB utilizza algoritmi di ricerca del cammino minimo per 
%      determinare tutti i percorsi più brevi tra coppie di nodi e quindi conteggia 
%      il numero di volte che ciascun nodo compare in questi percorsi.
%    - L'operazione viene eseguita per ogni nodo della rete, fornendo un quadro 
%      completo del suo ruolo come intermediario nei percorsi più brevi.

%    Limitazioni:
%    - La centralità di Betweenness può essere sensibile a modifiche nella 
%      struttura della rete, come l'aggiunta o rimozione di nodi o archi.
%    - Può non essere altrettanto informativa in reti con molti percorsi paralleli 
%      dove la funzione di Betweenness è distribuita tra più nodi.

% 4. Centralità Eigenvector(Eigenvector Centrality):
%    La centralità Eigenvector, valuta l'importanza 
%    di un nodo in una rete non solo in base al numero di collegamenti, ma anche 
%    in base alla qualità o all'importanza dei nodi a cui è collegato. Un nodo 
%    ha un'alta centralità di autovettore se è connesso a nodi che sono a loro 
%    volta importanti o ben connessi nella rete. 

%    Calcolo:
%    La centralità Eigenvector di un nodo in una rete è calcolata come il 
%    valore di un componente di un autovettore corrispondente all'autovalore 
%    massimo della matrice di adiacenza della rete. In termini matematici, se A 
%    è la matrice di adiacenza, allora la centralità di autovettore è data dall'equazione:
%    Ax = λx
%    dove A è la matrice di adiacenza, x è l'autovettore corrispondente all'autovalore 
%    massimo λ. Il valore di centralità per ogni nodo è dato dalla corrispondente 
%    componente nell'autovettore x.

%    Algoritmo MATLAB:
%    - MATLAB calcola la centralità Eigenvector utilizzando la funzione 'centrality'
%      con l'opzione 'eigenvector' su un oggetto grafo.
%    - Internamente, MATLAB utilizza metodi numerici per calcolare gli autovalori e 
%      gli autovettori della matrice di adiacenza(come il metodo dele potenze,QR,SVD). Il calcolo si concentra 
%      sull'identificazione dell'autovalore massimo e del suo corrispondente autovettore.
%    - L'autovettore risultante fornisce la misura della centralità di autovettore 
%      per ciascun nodo, con i valori maggiori che indicano un'alta centralità 
%      nell'ambito della rete.

%    Limitazioni:
%    - La centralità Eigenvector può essere influenzata da nodi che hanno molti 
%      collegamenti ma che non sono necessariamente centrali in termini di flusso 
%      di informazioni o risorse.
%    - Potrebbe non essere efficace in reti con molte componenti isolate o in reti 
%      dove la centralità è uniformemente distribuita.

% 5. Centralità PageRank:
%    La centralità PageRank è un'evoluzione della centralità Eigenvector
%    che introduce un fattore di casualità nel calcolo dell'importanza di un nodo 
%    in una rete. Inizialmente sviluppato da Google per classificare le pagine web, 
%    il PageRank assegna un punteggio di importanza a ciascun nodo basandosi non 
%    solo sul numero di collegamenti, ma anche sulla qualità e sull'importanza dei 
%    nodi collegati.

%    Calcolo:
%    Il calcolo del PageRank si basa sulla navigazione casuale di un 'surfer' in 
%    una rete. In ogni passo, il 'surfer' può scegliere di seguire un collegamento 
%    esistente con una certa probabilità (d, tipicamente 0.85) o saltare a un nodo 
%    casuale nella rete con probabilità 1-d. Il punteggio PageRank di un nodo è 
%    quindi la probabilità di trovarsi in quel nodo in un lungo termine di navigazione.
%    La formula del PageRank è:
%    PR(i) = (1-d)/N + d Σ (PR(j)/L(j)) per tutti i nodi j che puntano a i
%    dove N è il numero totale di nodi nella rete, PR(j) è il PageRank del nodo j, 
%    e L(j) è il numero di collegamenti uscenti dal nodo j.

%    Algoritmo MATLAB:
%    - MATLAB calcola il PageRank utilizzando la funzione 'centrality' con l'opzione 
%      'pagerank' su un oggetto grafo.
%    - L'algoritmo di PageRank implementato in MATLAB utilizza tecniche di algebra lineare 
%      e metodi iterativi per approssimare la distribuzione del PageRank su tutti 
%      i nodi della rete.
%    - Il processo iterativo continua finché la distribuzione di PageRank non converge 
%      a un valore stabile, indicando l'importanza relativa di ciascun nodo nella rete.

%    Limitazioni:
%    - In reti con molte componenti isolate o in reti non fortemente connesse, 
%      il PageRank può essere distorto, privilegiando nodi in componenti maggiormente connesse.
%    - La presenza di 'dead ends' (nodi senza collegamenti uscenti) o 'spider traps' 
%      (cicli di nodi che si riferiscono l'un l'altro) può influenzare il calcolo 
%      del PageRank e richiedere l'implementazione di strategie di correzione.

% 6. Centralità del Sottografo Esponenziale (Exponential Subgraph Centrality):
%    Questa forma di centralità si basa sull'idea che l'importanza di un nodo 
%    può essere determinata dal numero di sottografi esponenziali in cui è coinvolto.
%    Concettualmente, ciò richiederebbe il calcolo dell'esponenziale della matrice 
%    di adiacenza, un compito computazionalmente impegnativo.

%    In un contesto matematico più ampio, il calcolo dell'esponenziale di una matrice
%    può essere collegato agli integrali di Stieltjes. Questi integrali permettono 
%    di sommare in modo continuo le funzioni ponderate, fornendo una base teorica 
%    per considerare tutte le possibili connessioni (o percorsi) in una rete.

%    Tuttavia, il calcolo diretto di questi integrali è complesso, quindi si ricorre 
%    alle formule di quadratura per semplificarlo. Queste formule consentono di 
%    approssimare un integrale complesso attraverso somme finite, rendendo 
%    il calcolo più praticabile.

%    Per la centralità del sottografo esponenziale, si utilizza l'algoritmo di Lanczos, 
%    che opera sui sottospazi di Krylov. Questo metodo trasforma la matrice di 
%    adiacenza in una forma tridiagonale più semplice, mantenendo le informazioni 
%    essenziali sui percorsi nella rete. Il risultato è una matrice tridiagonale J, 
%    che è un'approssimazione semplificata ma efficace della matrice originale.

%    Il passo finale è il calcolo dell'esponenziale di J (expm(J)). Questo si avvicina 
%    all'esponenziale della matrice di adiacenza originale, catturando l'essenza 
%    della centralità del sottografo esponenziale in un modo computazionalmente 
%    efficiente.
%    Limitazioni:
%    - Computazionalmente onerosa, specialmente per reti di grandi dimensioni.
%    - Possibili problemi numerici e di stabilità con l'algoritmo di Lanczos.
%    - Potrebbe non essere informativa in reti con strutture irregolari o componenti isolate.


% 7. Centralità del Risolvente (Resolvent Centrality):
%    Il risolvente invece di fare somma righe come nell Katz Centrality prende l'elemento diagonale
%    Questo indice di centralità misura l'influenza di un nodo nella rete 
%    considerando non solo i collegamenti diretti, ma anche quelli indiretti 
%    di varie lunghezze. Utilizza la matrice risolvente, che è una forma 
%    di rappresentazione della rete che tiene conto dei cammini di lunghezza 
%    variabile tra i nodi.
%    
%    La centralità della risolvente di un nodo è influenzata dalla sua 
%    posizione nella rete e dalle connessioni dei nodi a cui è direttamente 
%    e indirettamente collegato. Un nodo con alta centralità della risolvente 
%    è posizionato in modo tale da poter influenzare efficacemente un gran 
%    numero di altri nodi nella rete, sia direttamente che indirettamente.
%    
%    Il calcolo si basa sull'inversione della matrice (I - αA), dove I è la 
%    matrice identità, A è la matrice di adiacenza della rete, e α è un 
%    parametro che determina l'importanza relativa dei percorsi a diverse 
%    distanze. Un valore di α piccolo dà più peso ai percorsi più lunghi, 
%    mentre un valore più alto enfatizza i collegamenti diretti.
%    Limitazioni:
%    - Richiede l'inversione della matrice (I - αA), che può essere dispendiosa e soggetta a errori numerici.
%    - Potrebbe non riflettere accuratamente l'importanza dei nodi in reti con strutture topologiche complesse.
%    - La scelta del parametro α è cruciale e può influenzare i risultati.


% 8. Centralità di Katz:
%    Fa somma righe. 
%    Generalizzazione della centralità eigenvector, che tiene conto anche dei percorsi indiretti.
%    Assegna un punteggio più alto ai nodi che sono connessi a molti nodi che sono a loro volta connessi.
%    La Centralità di Katz differisce dalla centralità Eigenvector per la sua capacità di considerare
%    anche percorsi di lunghezza maggiore, attribuendo un certo peso anche alle connessioni indirette.
%
%    Il calcolo della Centralità di Katz si basa sull'utilizzo di un parametro alpha (α), 
%    che modula l'importanza attribuita ai percorsi di diverse lunghezze. 
%    Un valore più basso di α dà maggiore importanza ai collegamenti indiretti.
%
% Due Metodi per il Calcolo della Centralità di Katz:
% - Metodo Lineare:
%   Questo metodo si basa sull'inversione della matrice (I - αA), dove I è la matrice identità 
%   e A è la matrice di adiacenza della rete. È più efficiente per reti di dimensioni ridotte 
%   o medie e meno sensibile a problemi numerici come l'overflow. Il risultato fornisce un punteggio 
%   di centralità che riflette l'influenza di un nodo tenendo conto sia dei collegamenti diretti 
%   che di quelli indiretti.
%
% - Metodo Basato su 'expm':
%   Utilizza la funzione matematica 'expm' per calcolare l'esponenziale della matrice di adiacenza shiftata (A - μI),
%   dove μ è generalmente scelto come l'autovalore massimo di A. Questo metodo è più adatto per reti di grandi 
%   dimensioni con una struttura complessa e molte connessioni, poiché cattura in modo più accurato l'influenza 
%   di percorsi lunghi e complessi. Tuttavia, può essere più lento e soggetto a problemi numerici come overflow
%   o imprecisioni nei calcoli degli esponenziali di matrici di grandi dimensioni.
%    Limitazioni:
%    - Computazionalmente intensiva per reti di grandi dimensioni, specialmente con il metodo 'expm'.
%    - Può essere influenzata da nodi ad alto grado che non sono centrali per il flusso di informazioni.
%    - La scelta di α può portare a risultati diversi, con rischio di sovrastima o sottovalutazione dei collegamenti indiretti.

% 9. Centralità di Reti Temporali (Temporal Network Centrality):
%    Valuta la centralità in reti che cambiano nel tempo, considerando 
%    la dinamica temporale nella valutazione della centralità dei nodi.
%    Si basa sull'uso di un tensore, che può esser visto come un'estensione multidimensionale 
%    di una matrice, per rappresentare le relazioni tra i nodi in diversi 
%    istanti temporali.

%    - Broadcast Centrality: 
%      Si calcola costruendo una matrice a blocchi (tensore) dove ogni blocco
%      rappresenta la matrice di adiacenza in un determinato istante temporale. 
%      La formula utilizzata è: B = [A1, 0, ..., 0; 0, A2, ..., 0; ...; 0, ..., 0, An],
%      dove Ai rappresenta la matrice di adiacenza all'istante i. La centralità 
%      broadcast di un nodo è data dalla somma dei valori nelle righe corrispondenti 
%      di questa matrice a blocchi, riflettendo la sua capacità di influenzare altri nodi.

%    - Receive Centrality: 
%      Simile al calcolo della broadcast centrality, ma in questo caso, si considera
%      la capacità di un nodo di ricevere informazioni. La matrice a blocchi è 
%      costruita invertendo l'ordine dei blocchi rispetto al caso broadcast.
%      La formula utilizzata è: B = [0, ..., 0, An; 0, ..., An-1, 0; ...; A1, 0, ..., 0].
%      La centralità receive è calcolata sommando i valori nelle colonne corrispondenti 
%      del nodo in questione nella matrice a blocchi.

%    In entrambi i casi, il calcolo richiede la costruzione di una matrice a blocchi
%    che rappresenta la rete in diversi momenti temporali, connessa in modo 
%    appropriato a seconda che si stia misurando la broadcast o la receive centrality.
%    - La scelta tra broadcast e receive dipende dall'obiettivo specifico dell'analisi:
%    se l'interesse è nel capire come le informazioni o influenze si diffondono dalla 
%    fonte, si utilizza la broadcast centrality; se l'interesse è nel capire come e 
%    quanto rapidamente un nodo riceve informazioni, si utilizza la receive centrality.    
%    Limitazioni:
%    - Gestione complessa e computazionalmente intensiva delle reti tramite tensori.
%    - Suscettibile a variazioni significative a seconda della scala temporale e dinamica della rete.
%    - Le misure di broadcast e receive centrality possono non catturare pienamente l'importanza dei nodi in momenti diversi.


% -- Comandi utili --
% Per ulteriori dettagli su ciascun indice, usare il comando `help` seguito dal nome dell'indice.
% Esempio: `help centrality`, help tnode_centralities.m per visualizzare questi commenti.
% close all,clear all,clc 

% Chiedi all'utente se desidera connettere il grafo
connectGraph = input('Vuoi aumentare la probabilità di avere una rete connessa? (si/no): ', 's');

% Genera il grafo complesso
n = randi([10, 200]);

if strcmpi(connectGraph, 'Si') || strcmpi(connectGraph, 'si') || strcmp(connectGraph, 's')
    % Scegliendo m_erdrey = ceil((n * log(n))/2), siamo vicino alla soglia critica dove 
    % la probabilità che il grafo sia connesso inizia a diventare significativa.
    m_erdrey = ceil((n * log(n))); % più probabile
else
    % Se l'utente non desidera un grafo connesso, si sceglie un valore di m inferiore,
    % mantenendo un numero di archi che non garantisce la connettività.
    m_erdrey = randi([10, 199]);
    %m_erdrey = 3; %ceil(n); 
end
    A = erdrey(n, m_erdrey);

% Crea il grafo
G = graph(A);

% Stampa del numero di nodi e archi
fprintf('Numero di nodi: %d\n', numnodes(G));
fprintf('Numero di archi: %d\n', numedges(G));

% Visualizzazione della rete
figure;
plot(G);
title('Visualizzazione della Rete Complessa');


% Generazione di M matrici di adiacenza temporali
M = 10;%randi([1, 10]);  % Numero di "fette" temporali, massimo 20
T = zeros(n, n, M);  % Inizializzazione del tensore T

while true
    choice = menu('Seleziona l''indice di centralità:', ...
        'Degree', 'Closeness', 'Betweenness', 'Eigenvector', 'PageRank', ...
        'Exp Subgraph', 'Resolvent', 'Katz', 'Temp Net', 'Color Plots', 'Esci');

    if choice == 11
        break; % Uscita dal ciclo se l'utente sceglie 'Esci'
    end

    % Richiesta del numero di nodi importanti
    m = input('Inserisci il numero di nodi importanti da identificare: ');
    while m > numnodes(G)
        fprintf('Il numero inserito è maggiore del numero totale di nodi (%d). Riprova.\n', numnodes(G));
        m = input('Inserisci il numero di nodi importanti da identificare: ');
    end

    if choice == 10
        % Chiamata alla funzione centr_colors
        centrality_colors(A, m); 
        continue; % Continua il ciclo per tornare al menu principale
    end
    % Calcolo dell'indice di centralità
    if choice <= 5
        centralityTypes = {'degree', 'closeness', 'betweenness', 'eigenvector', 'pagerank'};
        %tic
        centrality = centrality(G, centralityTypes{choice});
        %toc; % Termina il monitoraggio del tempo
        %fprintf('Tempo impiegato per %s centrality: %f secondi.\n', centralityTypes{choice}, toc);
       [sortedValues, sortedIndices] = sort(centrality, 'descend');
        topIndices = sortedIndices(1:m);
        topValues = sortedValues(1:m);
    elseif choice == 9
        close all;
        fprintf('Numero di "fette" temporali (%d).\n', M);
        n_temp = 10; % Dimensione fissa per le matrici di adiacenza temporali
        T_net = zeros(n_temp, n_temp, M); % Inizializzazione del tensore T_net
        
        for k = 1:M
            A_net = erdrey(n_temp, 9); % Generazione della matrice di dimensione n_temp
            A_net = A_net * A_net'; % Rendiamo la matrice simmetrica
            T_net(:,:,k) = A_net;
        end

         % Chiede all'utente di inserire un numero di nodi tra 2 e 9
        m_temp_net = input('Inserisci il numero di nodi importanti da identificare (tra 2 e 9): ');
        while m_temp_net < 2 || m_temp_net > 9
            fprintf('Numero non valido. Inserisci un numero tra 2 e 9.\n');
            m_temp_net = input('Inserisci il numero di nodi importanti da identificare (tra 2 e 9): ');
        end
        % Calcola e mostra i nodi di centralità
        centralityType = input('Scegli il tipo di centralità (broadcast/receive): ', 's');
        [topIndices, topValues] = temp_net(T_net, m, centralityType);
        
        % Visualizzazione multipla delle fette temporali
        numSubplotsPerRow = ceil(sqrt(M)); % Calcola il numero di subplot per riga
        figure; % Crea una nuova figura
        
        for k = 1:M
            A_net = T_net(:,:,k); % Estrai la k-esima matrice di adiacenza
            G_net = graph(A_net); % Crea il grafo per la k-esima fetta temporale
            
            subplot(numSubplotsPerRow, numSubplotsPerRow, k); % Crea un subplot per la k-esima fetta
            p = plot(G_net, 'Layout', 'force'); % Visualizza il grafo
            title(sprintf('Fetta Temporale %d', k)); % Aggiungi un titolo al subplot
            
            % Evidenzia i nodi importanti
            highlight(p, topIndices, 'NodeColor', 'r', 'MarkerSize', 10);
        end
        
        sgtitle('Visualizzazione Multipla delle Fette Temporali con Nodi Evidenziati'); % Titolo generale
        
        fprintf('Top %d nodi per l''indice di centralità scelto:\n', m);
        for i = 1:m
            fprintf('Indice %d: Centralità = %f\n', topIndices(i), topValues(i));
        end
        break;

    else
        switch choice
            case 6
                [topIndices, topValues] = exp_sub_centr(A, m);
            case 7
                [topIndices, topValues] = res_sub_centr(A, m);
            case 8
                fprintf('Scegli il metodo per calcolare Katz centrality:\n');
                fprintf('1 - Metodo lineare\n');
                fprintf('2 - Metodo basato su expm\n');
                methodChoice = input('Inserisci il numero corrispondente alla tua scelta: ');
                while methodChoice ~= 1 && methodChoice ~= 2
                    disp('Scelta non valida. Per favore scegli 1 per il metodo lineare o 2 per il metodo basato su expm.');
                    methodChoice = input('Inserisci il numero corrispondente alla tua scelta: ');
                end
                if methodChoice == 1
                    method = 'linear';
                else
                    method = 'exp';
                end
                [topIndices, topValues] = katz_centr(A, m, method);
        end
    end

    % Stampa dei risultati e preparazione dei dati per la visualizzazione grafica
    fprintf('Top %d nodi per l''indice di centralità scelto:\n', m);
    centralitiesForPlot = zeros(numnodes(G), 1); % Vettore per tutti i nodi inizializzato a zero
    
    for i = 1:m
        fprintf('Indice %d: Centralità = %f\n', topIndices(i), topValues(i));
        centralitiesForPlot(topIndices(i)) = topValues(i); % Imposta la centralità per i nodi selezionati
    end
    
    % Visualizzazione della rete con i nodi colorati in base alla loro centralità
    p = plot(G, 'MarkerSize', 6);
    p.NodeCData = centralitiesForPlot; % Assegna il vettore di centralità ai nodi nel grafico
    colormap jet; % Utilizza una colormap
    colorbar; % Mostra la colorbar per interpretare i valori di centralità
    title('Visualizzazione della Rete con Centralità dei Nodi');
end


% -- Funzioni di centralità --

% Funzione exp_sub_centr
function [i, val] = exp_sub_centr(A, m)
    n = size(A, 1); % Numero di nodi
    centralities = zeros(n, 1); % Vettore per memorizzare la centralità di ogni nodo
    kmax = min(n, 20); % Numero massimo di iterazioni
    U = zeros(n, kmax); % Preallocazione di U

    for node = 1:n
        u = zeros(n, 1);
        u(node) = 1; % Vettore iniziale
        J = zeros(n); % Matrice tridiagonale J
        tol = 1e-5; % Tolleranza
        g = 0; % Inizializzazione di g
        betaa = 0;% se non lo faccio da problemi
        
        for k = 1:kmax
            U(:, k) = u; % Aggiornamento di U
            alfa = u' * A * u;
            utilde = A * u - alfa * u - betaa * U(:, max(k-1, 1));
            betaa = norm(utilde);

            if betaa < tol
                break; % Uscita se beta è troppo piccolo
            end

            u = utilde / betaa;
            
            J(k, k) = alfa;
            if k < n
                J(k + 1, k) = betaa;
                J(k, k + 1) = betaa;
            end

            E = expm(J(1:k, 1:k));%Matrice esponenziale e non esponenziale dei coefficienti
            g = E(1, 1); % Centralità del nodo corrente

            if k > 1 && abs(g - E(1, 1)) < tol * abs(g)
                break;
            end
        end

        centralities(node) = g;
    end

    [sortedValues, sortedIndices] = sort(centralities, 'descend');
    i = sortedIndices(1:m);
    val = sortedValues(1:m);
end


function [i, val] = katz_centr(A, m, method)
%tic
    % Controllo del numero di argomenti in ingresso e impostazione del metodo di default.
    if nargin < 3
        method = 'linear'; % Metodo di default
    end
    
    % Dimensione della matrice di adiacenza.
    n = size(A, 1);
    % Calcolo del massimo autovalore in modulo per determinare alpha.
    rhoA = max(abs(eigs(A, 1)));
    % Impostazione di alpha per assicurare la convergenza (alpha < 1/rho(A)).
    alpha = 0.9 / rhoA;
    if strcmp(method, 'linear')
        % Metodo lineare: calcolo della centralità di Katz direttamente.
        assert(alpha < 1/rhoA, 'Alpha is too large and may cause the Katz centrality to diverge.');
        katzCentrality = (eye(n) - alpha * A) \ ones(n, 1);
    elseif strcmp(method, 'exp')
        % Metodo basato sulla funzione matriciale con uso di expm.
        
        % Lo shift della matrice è una tecnica numerica per aumentare la stabilità
        % del calcolo dell'esponenziale di una matrice, particolarmente utile
        % per evitare l'overflow in matrici con autovalori di grande modulo.
        mu = rhoA; % Autovalore massimo per lo shift.
        A_shifted = A - mu * eye(n); % Applicazione dello shift.
        
        % Uso della funzione expm per calcolare l'esponenziale della matrice shiftata.
        % expm è basata su algoritmi che sfruttano i sottospazi di Krylov
        % per una computazione efficiente e accurata dell'esponenziale di matrice.
        try
            katzCentrality = expm(A_shifted) * ones(n, 1);
        catch
            error('Overflow encountered in matrix exponential');
        end
        
        % Dopo aver calcolato l'esponenziale di A_shifted, è necessario
        % "disfare" lo shift moltiplicando per e^mu.
        katzCentrality = katzCentrality * exp(mu);
        
        % Normalizzazione per confronto equo dei valori di centralità.
        katzCentrality = katzCentrality / max(katzCentrality);
    else
        % Gestione dell'errore se il metodo fornito non è riconosciuto.
        error('Metodo non riconosciuto. Scegli "linear" o "exp".');
    end
    
    % Ordinamento dei nodi in base al loro valore di centralità e selezione dei top m.
    [sortedValues, sortedIndices] = sort(katzCentrality, 'descend');
    i = sortedIndices(1:m);
    val = sortedValues(1:m);
    %toc
%fprintf('Tempo impiegato per katz_centr: %f secondi.\n', toc);
end

function [i, val] = res_sub_centr(A, m)
%tic
    n = size(A, 1);
    rhoA = max(abs(eigs(A, 1)));
    alpha = 0.9 / rhoA;
    resolventCentrality = zeros(n, 1);
    for j = 1:n
        ei = zeros(n, 1);
        ei(j) = 1;
        y = (eye(n) - alpha * A) \ ei;
        resolventCentrality(j) = y(j);
    end
    [sortedValues, sortedIndices] = sort(resolventCentrality, 'descend');
    i = sortedIndices(1:m);
    val = sortedValues(1:m);
%toc
%fprintf('Tempo impiegato per res_sub_centr: %f secondi.\n', toc);
 end

function [i, val] = temp_net(T, m, centralityType)
tic
    M = size(T, 3);  % Numero di fette temporali
    n = size(T, 1);  % Dimensione di ciascuna rete
    B = zeros(M*n);  % Inizializzazione della matrice a blocchi

    % Calcolo del più grande autovalore in valore assoluto per ogni "fetta" temporale
    rho = max(arrayfun(@(k) abs(eigs(T(:,:,k), 1, 'largestabs')), 1:M));
    alpha = 0.9 / rho;

    % Riempimento della matrice a blocchi
    for k = 1:M
        idx = (k-1)*n + (1:n);  % Indici per la k-esima fetta temporale

        B(idx, idx) = alpha * T(:,:,k);

        if k < M
            if strcmp(centralityType, 'broadcast')
                B(idx, idx+n) = eye(n);
            elseif strcmp(centralityType, 'receive')
                B(idx+n, idx) = eye(n);
            end
        end
    end

    if strcmp(centralityType, 'broadcast')
        E = [zeros((M-1)*n, 1); ones(n, 1)];  % Per la centralità broadcast
    elseif strcmp(centralityType, 'receive')
        E = [ones(n, 1); zeros((M-1)*n, 1)];  % Per la centralità receive
    end

    xbig = (eye(M*n) - B) \ E;
    xbig = xbig(1:n);  % Consideriamo solo la prima "fetta" temporale

    [sortedValues, sortedIndices] = sort(xbig, 'descend');
    i = sortedIndices(1:m);
    val = sortedValues(1:m);
    toc
fprintf('Tempo impiegato per temp_net: %f secondi.\n', toc);
end

function centrality_colors(A, m)
    % Creazione del grafo da una matrice di adiacenza A
    G = graph(A, 'omitselfloops');
    
    % Inizia la misura del tempo totale
    total_start_time = tic;

    % Calcolo delle diverse centralità
    centrDegree = centrality(G, 'degree');
    centrCloseness = centrality(G, 'closeness');
    centrBetweenness = centrality(G, 'betweenness');
    centrEigenvector = centrality(G, 'eigenvector');
    centrPageRank = centrality(G, 'pagerank');

    % Inizializza una cell array per memorizzare i risultati
    results = cell(5, 4); % Aggiunta di una colonna per i valori dei nodi migliori
    
    % Funzione per plottare il grafo e evidenziare i top m nodi
    function plotAndHighlight(graph, centralityValues, subplotIndex, numTopNodes, titleText)
        subplot_start_time = tic; % Inizia la misurazione del tempo per il subplot
        subplot(2, 3, subplotIndex);
        p = plot(graph, 'MarkerSize', 4);
        p.NodeCData = centralityValues;
        colormap turbo;
        colorbar;
        title(titleText);

        % Ordina e seleziona i top m nodi
        [~, sortedIndices] = sort(centralityValues, 'descend');
        topIndices = sortedIndices(1:numTopNodes);

        % Evidenzia i top m nodi con un marker più grande
        highlight(p, topIndices, 'MarkerSize', 10);

        % Memorizza i risultati nella cell array
        topValues = centralityValues(topIndices);
        results{subplotIndex, 1} = titleText;
        results{subplotIndex, 2} = topIndices;
        results{subplotIndex, 3} = topValues;
        results{subplotIndex, 4} = toc(subplot_start_time); % Tempo impiegato per questa centralità
        
        % Stampa il tempo per questa centralità sotto il subplot
        subplot(2, 3, subplotIndex);
        text(0.5, -0.15, sprintf('Time: %.4f sec', results{subplotIndex, 4}), 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end

    % Plotta e evidenzia i top m nodi per ciascuna centralità
    plotAndHighlight(G, centrDegree, 1, m, 'Degree Centrality');
    plotAndHighlight(G, centrCloseness, 2, m, 'Closeness Centrality');
    plotAndHighlight(G, centrBetweenness, 3, m, 'Betweenness Centrality');
    plotAndHighlight(G, centrEigenvector, 4, m, 'Eigenvector Centrality');
    plotAndHighlight(G, centrPageRank, 5, m, 'PageRank Centrality');

    % Visualizzazione del tempo totale e stampa della tabella
    total_end_time = toc(total_start_time);
    fprintf('\nTop %d Nodi per ogni Centralità e Tempo Impiegato:\n', m);
    fprintf('-------------------------------------------------------------\n');
    fprintf('%-25s %-15s %-15s %-15s\n', 'Centralità', 'Nodi Migliori', 'Tempo (sec)', 'Tempo per Calcolo');
    fprintf('-------------------------------------------------------------\n');
    for i = 1:5
        centralityName = results{i, 1};
        topNodes = results{i, 2};
        elapsedTime = results{i, 4};
        computationTime = toc(total_start_time) - elapsedTime;
        topNodesStr = strjoin(arrayfun(@num2str, topNodes, 'UniformOutput', false), ', ');
        fprintf('%-25s %-15s %-15.4f %-15.4f\n', centralityName, topNodesStr, elapsedTime, computationTime);
    end
    
    % Stampa dei valori delle centralità
    fprintf('\nValori di Centralità per i Nodi Migliori:\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('%-25s %-15s\n', 'Centralità', 'Valori');
    fprintf('-------------------------------------------------------------\n');
    for i = 1:5
        centralityName = results{i, 1};
        topValues = results{i, 3};
        topValuesStr = strjoin(arrayfun(@(x) sprintf('%.4f', x), topValues, 'UniformOutput', false), ', ');
        fprintf('%-25s %-15s\n', centralityName, topValuesStr);
    end
    
    % Stampa del tempo totale impiegato
    fprintf('\nTempo totale: %.4f sec\n', total_end_time);
end



