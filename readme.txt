Plik dane/ciala.txt zawiera informacje o cz³onach uk³adu.
W pierwszej linijce - iloœæ cz³onów. W kolejnych - informacje o cz³onach: x-owa, y-owa wspó³rzêdna œrodka masy cz³onu oraz k¹t fi.

Plik dane/pary.txt zawiera informacje o parach kinematycznych w uk³adzie.
W pierwszej linijce - iloœæ par obrotowych. W nastêpnych - definicje par: indeksy cz³onów tworz¹cych parê i, j (kolejnoœæ jak w dane/ciala.txt); lokalizacja pary obrotowej w uk³adzie globalnym.
Nastêpnie iloœæ par postêpowych i kolejno ich definicje: indeksy cz³onów i,j; punkt definiuj¹cy parê w cz³onie i; punkt definiuj¹cy parê w cz³onie j. Wspó³rzêdne w uk³adzie globalnym.

Plik dane/wymuszenia.txt zawiera informacje o wymuszeniach.
W pierwszej linijce - iloœæ wymuszeñ w parach obrotowych. W kolejnych: indeks pary obrotowej (kolejnoœæ jak w dane/pary.txt) i indeks wymuszenia.
Nastêpnie - iloœæ wymuszeñ w parach postêpowych. W kolejnych: indeks pary postêpowej (kolejnoœæ jak w dane/pary.txt) i indeks wymuszenia.
Definicje wymuszeñ znajduj¹ siê w plikach: Wymuszenie.m, DWymuszenie.m oraz DDWymuszenie.m. Indeksy wymuszeñ w instrukcji switch - case jak w pliku dane/wymuszenia.txt.

Plik doexportu.txt zawiera informacje o punktach których po³o¿enia, prêdkoœci i przyspieszenia maj¹ byæ zwracane. W pierwszej linijce iloœæ punktów, a nastêpnie dla ka¿dego z nich: indeks cz³onu do którego nale¿y oraz jego wspó³rzêdne w uk³adzie globalnym.