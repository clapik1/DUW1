Plik dane/ciala.txt zawiera informacje o cz�onach uk�adu.
W pierwszej linijce - ilo�� cz�on�w. W kolejnych - informacje o cz�onach: x-owa, y-owa wsp�rz�dna �rodka masy cz�onu oraz k�t fi.

Plik dane/pary.txt zawiera informacje o parach kinematycznych w uk�adzie.
W pierwszej linijce - ilo�� par obrotowych. W nast�pnych - definicje par: indeksy cz�on�w tworz�cych par� i, j (kolejno�� jak w dane/ciala.txt); lokalizacja pary obrotowej w uk�adzie globalnym.
Nast�pnie ilo�� par post�powych i kolejno ich definicje: indeksy cz�on�w i,j; punkt definiuj�cy par� w cz�onie i; punkt definiuj�cy par� w cz�onie j. Wsp�rz�dne w uk�adzie globalnym.

Plik dane/wymuszenia.txt zawiera informacje o wymuszeniach.
W pierwszej linijce - ilo�� wymusze� w parach obrotowych. W kolejnych: indeks pary obrotowej (kolejno�� jak w dane/pary.txt) i indeks wymuszenia.
Nast�pnie - ilo�� wymusze� w parach post�powych. W kolejnych: indeks pary post�powej (kolejno�� jak w dane/pary.txt) i indeks wymuszenia.
Definicje wymusze� znajduj� si� w plikach: Wymuszenie.m, DWymuszenie.m oraz DDWymuszenie.m. Indeksy wymusze� w instrukcji switch - case jak w pliku dane/wymuszenia.txt.

Plik doexportu.txt zawiera informacje o punktach kt�rych po�o�enia, pr�dko�ci i przyspieszenia maj� by� zwracane. W pierwszej linijce ilo�� punkt�w, a nast�pnie dla ka�dego z nich: indeks cz�onu do kt�rego nale�y oraz jego wsp�rz�dne w uk�adzie globalnym.