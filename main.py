import time
import copy
import pygame
import sys

ADANCIME_MAX = 4




class Joc:
    """
    Clasa care defineste jocul. Se va schimba de la un joc la altul.
    """
    celuleGrid = None
    JMIN = None
    JMAX = None
    GOL = '#'
    NR_LINII = 10
    NR_COLOANE = 10
    scor_x = None
    scor_0 = None

    def __init__(self, matr=None, NR_LINII=None, NR_COLOANE=None):
        # creez proprietatea ultima_mutare # (l,c)
        self.ultima_mutare = None

        if matr:
            # e data tabla, deci suntem in timpul jocului
            self.matr = matr
        else:
            # nu e data tabla deci suntem la initializare
            self.matr = [[self.__class__.GOL] * NR_COLOANE for i in range(NR_LINII)]

            if NR_LINII is not None:
                self.__class__.NR_LINII = NR_LINII
            if NR_COLOANE is not None:
                self.__class__.NR_COLOANE = NR_COLOANE



    def deseneaza_grid(self, coloana_marcaj=None):  # tabla de exemplu este ["#","x","#","0",......]

        for ind in range(self.__class__.NR_COLOANE * self.__class__.NR_LINII):
            linie = ind // self.__class__.NR_COLOANE  # // inseamna div
            coloana = ind % self.__class__.NR_COLOANE

            if coloana == coloana_marcaj:
                # daca am o patratica selectata, o desenez cu rosu
                culoare = (255, 255, 0)
            else:
                # altfel o desenez cu alb
                culoare = (255, 255, 255)
            pygame.draw.rect(self.__class__.display, culoare, self.__class__.celuleGrid[ind])  # alb = (255,255,255)
            if self.matr[linie][coloana] == 'x':
                self.__class__.display.blit(self.__class__.x_img, (
                coloana * (self.__class__.dim_celula + 1), linie * (self.__class__.dim_celula + 1)))
            elif self.matr[linie][coloana] == '0':
                self.__class__.display.blit(self.__class__.zero_img, (
                coloana * (self.__class__.dim_celula + 1), linie * (self.__class__.dim_celula + 1)))
        pygame.display.flip()
        #pygame.display.update()

    @classmethod
    def jucator_opus(cls, jucator):
        return cls.JMAX if jucator == cls.JMIN else cls.JMIN

    @classmethod
    def initializeaza(cls, display, NR_LINII=10, NR_COLOANE=10, dim_celula=15):
        cls.display = display
        cls.dim_celula = dim_celula
        cls.x_img = pygame.image.load('ics.png')
        cls.x_img = pygame.transform.scale(cls.x_img, (dim_celula, dim_celula))
        cls.zero_img = pygame.image.load('zero.png')
        cls.zero_img = pygame.transform.scale(cls.zero_img, (dim_celula, dim_celula))
        cls.celuleGrid = []  # este lista cu patratelele din grid
        for linie in range(NR_LINII):
            for coloana in range(NR_COLOANE):
                patr = pygame.Rect(coloana * (dim_celula + 1), linie * (dim_celula + 1), dim_celula, dim_celula)
                cls.celuleGrid.append(patr)

    def parcurgere(self, directie):
        um = self.ultima_mutare  # (l,c)
        culoare = self.matr[um[0]][um[1]]
        nr_mutari = 0
        while True:
            um = (um[0] + directie[0], um[1] + directie[1])
            if not 0 <= um[0] < self.__class__.NR_LINII or not 0 <= um[1] < self.__class__.NR_COLOANE:
                break
            if not self.matr[um[0]][um[1]] == culoare:
                break
            nr_mutari += 1
        return nr_mutari


    def final(self): #punctul 7 din barem
        if(self.__class__.GOL not in self.matr or self.matr.count(self.__class__.GOL) ==1) :
            print("Nu se mai pot face mutari, jocul s-a terminat!")
            return True

        for i in range(len(self.matr)):
            for j in range(len(self.matr)):
                linie_i = i//10
                coloana_i = i%10
                linie = j//10
                coloana = j%10
                if(self.matr[i] == self.__class__.GOL and self.matr[j] == self.__class__.GOL): #am spatii libere ca sa pot face mutare
                    if(self.mutare_valida(linie_i, coloana_i) and (self.pot_completa_piesa(linie_i, coloana_i,linie,coloana))):
                        return False #mergem mai departe, nu se incheie
                        #verific in plus daca e o mutare valida si daca pot pune si a doua piesa

        return True


    def calculeaza_scor(self):
        Joc.scor_x=0
        Joc.scor_0=0

        for i in range(10):
            for j in range(10):
                if((self.matr[i*10 + j]) != self.__class__.GOL  and ((self.matr[i*10 + j]) == self.matr[(i-1)*10 + j+1]  )
                    and self.matr[(i-1)*10 + j+1] == self.matr[(i-2)*10 + j+2]): #diag dreapta sus
                    if(self.matr[i*10 + j] == "x"):
                        Joc.scor_x += 1

                    if(self.matr[i*10 + j] == "0"):
                        Joc.scor_0 += 1

                elif( (self.matr[i*10 + j])!= self.__class__.GOL  and self.matr[i*10 + j] == self.matr[(i+1)*10 + j+1]
                    and self.matr[(i+1)*10 + j+1] == self.matr[(i+2)*10 + j+2]):    #diag dreapta jos
                    if (self.matr[i * 10 + j] == "x"):
                        Joc.scor_x += 1

                    if (self.matr[i * 10 + j] == "0"):
                        Joc.scor_0 += 1


                elif( (self.matr[i*10 + j]) != self.__class__.GOL and self.matr[i*10 + j]==self.matr[(i-1)*10 + j-1]
                    and self.matr[(i-1)*10 + j-1] == self.matr[(i-2)*10 + j-2]):  #stanga sus
                    if (self.matr[i * 10 + j] == "x"):
                        Joc.scor_x += 1

                    if (self.matr[i * 10 + j] == "0"):
                        Joc.scor_0 += 1
                elif( (self.matr[i*10 + j]) != self.__class__.GOL and self.matr[i*10 + j]==self.matr[(i+1)*10 + j-1]
                    and self.matr[(i+1)*10 + j-1] == self.matr[(i+2)*10+ j-1]): #stanga jos
                    if (self.matr[i * 10 + j] == "x"):
                        Joc.scor_x += 1

                    if (self.matr[i * 10 + j] == "0"):
                        Joc.scor_0 += 1






    def estimeaza_scor(self, adancime):     #punctul 8
        final = self.final()
        self.calculeaza_scor()

        if final:
            if(Joc.scor_x == Joc.scor_0) :
                castigator = "remiza"

            elif Joc.scor_x > Joc.scor_0:
                castigator = "x"

            else:
                castigator = "0"

            if(castigator == Joc.JMAX):
                return (90+adancime)

            elif (castigator == Joc.JMIN):
                return (-90-adancime)

            else:
                return 0 #avem remiza

        elif self.__class__.JMAX == "x" :
            return Joc.scor_x - Joc.scor_0
        else:
            return Joc.scor_0 - Joc.scor_x


    def vecini_linie(self, linie, coloana):
        if( (self.matr[linie]*10 +coloana+1)=="x" or (self.matr[linie]*10 +coloana+1)=="0"
            or (self.matr[linie]*10 +coloana-1)=="x" or (self.matr[linie]*10 +coloana-1)=="0"):
            return True
        return False

    def vecini_coloane(self,linie,coloana):
        if((self.matr[linie-1]*10 + coloana) == "x" or (self.matr[linie-1]*10 +coloana) == "0" or
                (self.matr[linie+1]*10 +coloana) == "x" or (self.matr[linie+1]*10 + coloana =="0")
        ):
            return True
        return False

    def vecini_diagonale(self,linie, coloana):
        if((self.matr[linie-1]*10 + coloana -1) == "x" or (self.matr[linie-1]*10+ coloana-1) == "0" or
                (self.matr[linie+1]*10 + coloana+1)=="x" or (self.matr[linie+1]*10+coloana+1)=="0" or
                (self.matr[linie+1]*10 + coloana-1)=="x" or (self.matr[linie+1]*10+coloana-1)=="0" or
                (self.matr[linie-1]*10 +coloana+1)=="x" or (self.matr[linie-1]*10+coloana+1)=="0"):
            return True
        return False

    def obtine_lista_vecini(self, linie, coloana):
        lista_vecini=[]

        if((linie*10 +coloana)%10 != 0):       #daca nu am plasat pe prima coloana, caut vecin in stanga
            vec_stanga = (linie*10 + coloana -1)
            lista_vecini.append(self.matr[vec_stanga])

        if((linie*10 + coloana)%10 !=0 and (linie*10+coloana)>9): #daca nu am pe prima linie sau prima col caut stanga sus
            vec_stanga_sus = ((linie-1)*10 + coloana-1)
            lista_vecini.append(self.matr[vec_stanga_sus])

        if((linie*10 + coloana)>9): #daca nu am pe prima linie, caut deasupra
            vec_deasupra = ((linie-1)*10 +coloana)
            lista_vecini.append(self.matr[vec_deasupra])

        if((linie*10 +coloana) > 9 and (linie*10+coloana)%10 != 9): #daca nu e pe prima linie si ultima coloana
            vec_dreapta_sus = ((linie-1)*10 + coloana +1)
            lista_vecini.append(self.matr[vec_dreapta_sus])

        if((linie*10 + coloana)%10 != 9): #daca nu e pe ultima coloana
            vec_dreapta = (linie*10 +coloana + 1)
            lista_vecini.append(self.matr[vec_dreapta])

        if((linie*10 + coloana)<90 and (linie*10 +coloana)%10 != 9):
            vec_dreapta_jos = ((linie+1)*10 + coloana+1)
            lista_vecini.append(self.matr[vec_dreapta_jos])

        if((linie*10+coloana)<90):
            vec_dedesubt = ((linie+1)*10 + coloana)
            lista_vecini.append(self.matr[vec_dedesubt])

        if((linie*10+coloana)<90 and (linie*10+coloana)%10 !=0):
            vec_stanga_jos = ((linie+1)*10+coloana-1)
            lista_vecini.append(vec_stanga_jos)

        return lista_vecini





    def mutare_valida(self,linie,coloana):
        lista_vecini = self.obtine_lista_vecini(linie,coloana)
        if 'x' in lista_vecini and '0' in lista_vecini :
            return True
        return False

    def pot_completa_piesa(self,linie_jumate_piesa, coloana_jumate_piesa, linie, coloana):
        if((linie_jumate_piesa == linie  and coloana_jumate_piesa == coloana+1)  #dreapta
            or (linie_jumate_piesa == linie-1 and coloana_jumate_piesa == coloana) #jos
            or (linie_jumate_piesa == linie and coloana_jumate_piesa == coloana-1) #stanga
            or (linie_jumate_piesa == linie and coloana_jumate_piesa == coloana+1)): #sus
            return True
        return False





    def mutari(self, jucator):   #punctul 5 barem
        #lista de mutari valide ale calculatorului
        l_mutari = []
        for i in range(len(self.matr)):
            for j in range(len(self.matr)):
                if(self.matr[i] == self.__class__.GOL and self.matr[j] == self.__class__GOL): #daca avem pozitii goale
                    linie_i = i//10
                    coloana_i = i%10
                    linie_j = j//10
                    coloana_j = j%10
                    if(self.matr[i]==Joc.GOL and self.matr[j]==Joc.GOL) \
                        and self.mutare_valida(linie_i, coloana_i) \
                        and self.pot_completa_piesa(linie_i, coloana_i, linie_j, coloana_j):
                        matr_tabla_noua = list(self.matr)
                        matr_tabla_noua[i] = jucator
                        matr_tabla_noua[j] = jucator
                        l_mutari.append(matr_tabla_noua)
        return l_mutari








    def sirAfisare(self):
        sir = "  |"
        sir += " ".join([str(i) for i in range(self.NR_COLOANE)]) + "\n"
        sir += "-" * (self.NR_COLOANE + 1) * 2 + "\n"
        sir += "\n".join([str(i) + " |" + " ".join([str(x) for x in self.matr[i]]) for i in range(len(self.matr))])
        return sir

    def __str__(self):
        return self.sirAfisare()

    def __repr__(self):
        return self.sirAfisare()




class Stare:
    """
    Clasa folosita de algoritmii minimax si alpha-beta
    Are ca proprietate tabla de joc
    Functioneaza cu conditia ca in cadrul clasei Joc sa fie definiti JMIN si JMAX (cei doi jucatori posibili)
    De asemenea cere ca in clasa Joc sa fie definita si o metoda numita mutari() care ofera lista cu configuratiile posibile in urma mutarii unui jucator
    """

    def __init__(self, tabla_joc, j_curent, adancime, parinte=None, scor=None):
        self.tabla_joc = tabla_joc
        self.j_curent = j_curent

        # adancimea in arborele de stari
        self.adancime = adancime

        # scorul starii (daca e finala) sau al celei mai bune stari-fiice (pentru jucatorul curent)
        self.scor = scor

        # lista de mutari posibile din starea curenta
        self.mutari_posibile = []

        # cea mai buna mutare din lista de mutari posibile pentru jucatorul curent
        self.stare_aleasa = None

    def mutari(self):   #aici avem mutarile posibile pe care le poate face jucatorul
        l_mutari = self.tabla_joc.mutari(self.j_curent)
        juc_opus = Joc.jucator_opus(self.j_curent)
        l_stari_mutari = [Stare(mutare, juc_opus, self.adancime - 1, parinte=self) for mutare in l_mutari]

        return l_stari_mutari

    def __str__(self):
        sir = str(self.tabla_joc) + "(Juc curent:" + self.j_curent + ")\n"
        return sir

    def __repr__(self):
        sir = str(self.tabla_joc) + "(Juc curent:" + self.j_curent + ")\n"
        return sir


""" Algoritmul MinMax """


def min_max(stare):
    if stare.adancime == 0 or stare.tabla_joc.final():
        stare.scor = stare.tabla_joc.estimeaza_scor(stare.adancime)
        return stare

    # calculez toate mutarile posibile din starea curenta
    stare.mutari_posibile = stare.mutari()

    # aplic algoritmul minimax pe toate mutarile posibile (calculand astfel subarborii lor)
    mutari_scor = [min_max(mutare) for mutare in stare.mutari_posibile]

    if stare.j_curent == Joc.JMAX:
        # daca jucatorul e JMAX aleg starea-fiica cu scorul maxim
        stare.stare_aleasa = max(mutari_scor, key=lambda x: x.scor)
    else:
        # daca jucatorul e JMIN aleg starea-fiica cu scorul minim
        stare.stare_aleasa = min(mutari_scor, key=lambda x: x.scor)
    stare.scor = stare.stare_aleasa.scor
    return stare


def alpha_beta(alpha, beta, stare):
    if stare.adancime == 0 or stare.tabla_joc.final():
        stare.scor = stare.tabla_joc.estimeaza_scor(stare.adancime)
        return stare

    if alpha > beta:
        return stare  # este intr-un interval invalid deci nu o mai procesez

    stare.mutari_posibile = stare.mutari()

    if stare.j_curent == Joc.JMAX:
        scor_curent = float('-inf')

        for mutare in stare.mutari_posibile:
            # calculeaza scorul
            stare_noua = alpha_beta(alpha, beta, mutare)

            if (scor_curent < stare_noua.scor):
                stare.stare_aleasa = stare_noua
                scor_curent = stare_noua.scor
            if (alpha < stare_noua.scor):
                alpha = stare_noua.scor
                if alpha >= beta:
                    break

    elif stare.j_curent == Joc.JMIN:
        scor_curent = float('inf')

        for mutare in stare.mutari_posibile:

            stare_noua = alpha_beta(alpha, beta, mutare)

            if (scor_curent > stare_noua.scor):
                stare.stare_aleasa = stare_noua
                scor_curent = stare_noua.scor

            if (beta > stare_noua.scor):
                beta = stare_noua.scor
                if alpha >= beta:
                    break
    stare.scor = stare.stare_aleasa.scor

    return stare


def afis_daca_final(stare_curenta):     #punctul 1
    final = stare_curenta.tabla_joc.final()
    if (final):
        if (final == "remiza"):
            print("Remiza!")
        else:
            print("A castigat " + final)

        return True

    return False


class Buton:
    def __init__(self, display=None, left=0, top=0, w=0, h=0, culoareFundal=(53, 80, 115),
                 culoareFundalSel=(89, 134, 194), text="", font="arial", fontDimensiune=16, culoareText=(255, 255, 255),
                 valoare=""):
        self.display = display
        self.culoareFundal = culoareFundal
        self.culoareFundalSel = culoareFundalSel
        self.text = text
        self.font = font
        self.w = w
        self.h = h
        self.selectat = False
        self.fontDimensiune = fontDimensiune
        self.culoareText = culoareText
        # creez obiectul font
        fontObj = pygame.font.SysFont(self.font, self.fontDimensiune)
        self.textRandat = fontObj.render(self.text, True, self.culoareText)
        self.dreptunghi = pygame.Rect(left, top, w, h)
        # aici centram textul
        self.dreptunghiText = self.textRandat.get_rect(center=self.dreptunghi.center)
        self.valoare = valoare

    def selecteaza(self, sel):
        self.selectat = sel
        self.deseneaza()

    def selecteazaDupacoord(self, coord):
        if self.dreptunghi.collidepoint(coord):
            self.selecteaza(True)
            return True
        return False

    def updateDreptunghi(self):
        self.dreptunghi.left = self.left
        self.dreptunghi.top = self.top
        self.dreptunghiText = self.textRandat.get_rect(center=self.dreptunghi.center)

    def deseneaza(self):
        culoareF = self.culoareFundalSel if self.selectat else self.culoareFundal
        pygame.draw.rect(self.display, culoareF, self.dreptunghi)
        self.display.blit(self.textRandat, self.dreptunghiText)


class GrupButoane:
    def __init__(self, listaButoane=[], indiceSelectat=0, spatiuButoane=10, left=0, top=0):
        self.listaButoane = listaButoane
        self.indiceSelectat = indiceSelectat
        self.listaButoane[self.indiceSelectat].selectat = True
        self.top = top
        self.left = left
        leftCurent = self.left
        for b in self.listaButoane:
            b.top = self.top
            b.left = leftCurent
            b.updateDreptunghi()
            leftCurent += (spatiuButoane + b.w)

    def selecteazaDupacoord(self, coord):
        for ib, b in enumerate(self.listaButoane):
            if b.selecteazaDupacoord(coord):
                self.listaButoane[self.indiceSelectat].selecteaza(False)
                self.indiceSelectat = ib
                return True
        return False

    def deseneaza(self):
        # atentie, nu face wrap
        for b in self.listaButoane:
            b.deseneaza()

    def getValoare(self):
        return self.listaButoane[self.indiceSelectat].valoare


############# ecran initial ########################
def deseneaza_alegeri(display, tabla_curenta):
    btn_alg = GrupButoane(
        top=30,
        left=30,
        listaButoane=[
            Buton(display=display, w=80, h=30, text="minimax", valoare="minimax"),
            Buton(display=display, w=80, h=30, text="alphabeta", valoare="alphabeta")
        ],
        indiceSelectat=1)
    btn_juc = GrupButoane(
        top=100,
        left=30,
        listaButoane=[
            Buton(display=display, w=35, h=30, text="x", valoare="x"),
            Buton(display=display, w=35, h=30, text="zero", valoare="0")
        ],
        indiceSelectat=0)

    btn_dificultate = GrupButoane(
        top=150,
        left=30,
        listaButoane=[
            Buton(display=display, w=80, h=30, text="incepator", valoare="incepator"), #punctul 2
            Buton(display=display, w=80, h=30, text="normal", valoare="medium"),
            Buton(display=display, w=80, h=30, text="avansat", valoare="avansat")
        ], indiceSelectat=0
    )
    ok = Buton(display=display, top=200, left=30, w=40, h=30, text="ok", culoareFundal=(155, 0, 55))
    btn_alg.deseneaza()
    btn_juc.deseneaza()
    btn_dificultate.deseneaza()
    ok.deseneaza()
    while True:
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            elif ev.type == pygame.MOUSEBUTTONDOWN:
                pos = pygame.mouse.get_pos()
                if not btn_alg.selecteazaDupacoord(pos):
                    if not btn_juc.selecteazaDupacoord(pos):
                        if not btn_dificultate.selecteazaDupacoord(pos):
                            if ok.selecteazaDupacoord(pos):
                                display.fill((0, 0, 0))  # stergere ecran
                                tabla_curenta.deseneaza_grid()
                                return btn_juc.getValoare(), btn_alg.getValoare(), btn_dificultate.getValoare()
        pygame.display.update()


def main():
    # setari interf grafica
    global jumate_piesa_plasata
    pygame.init()
    pygame.display.set_caption("RADU_GEORGE_DANIEL_311")
    # dimensiunea ferestrei in pixeli
    nl = 10
    nc = 10
    w = 50
    ecran = pygame.display.set_mode(size=(nc * (w + 1) - 1, nl * (w + 1) - 1))  # N *w+ N-1= N*(w+1)-1
    Joc.initializeaza(ecran, NR_LINII=nl, NR_COLOANE=nc, dim_celula=w)

    # initializare tabla
    tabla_curenta = Joc(NR_LINII=10, NR_COLOANE=10);
    #preluare butoane
    Joc.JMIN, tip_algoritm, dif  = deseneaza_alegeri(ecran, tabla_curenta)


    #Simbolul calculatorului
    Joc.JMAX = '0' if Joc.JMIN == 'x' else 'x'
    print("Simbol utilizator: ",Joc.JMIN)
    print("Algoritm: ", tip_algoritm)
    print("Dificultate: ", dif)

    print("Tabla initiala")
    print(str(tabla_curenta))

    #setare dificultate
    ADANCIME_MAX = 1 if dif == "incepator" else 2 if dif == "medium" else 3 #punctul 2

    # creare stare initiala
    stare_curenta = Stare(tabla_curenta, 'x', ADANCIME_MAX)

    tabla_curenta.deseneaza_grid()

    while True:

        if (stare_curenta.j_curent == Joc.JMIN):

            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    # iesim din program
                    pygame.quit()
                    sys.exit()
                if event.type == pygame.MOUSEMOTION:

                    pos = pygame.mouse.get_pos()  # coordonatele cursorului
                    for np in range(len(Joc.celuleGrid)):
                        if Joc.celuleGrid[np].collidepoint(pos):
                            stare_curenta.tabla_joc.deseneaza_grid(coloana_marcaj=np % Joc.NR_COLOANE)
                            break

                elif event.type == pygame.MOUSEBUTTONDOWN:

                    pos = pygame.mouse.get_pos()  # coordonatele cursorului la momentul clickului

                    for np in range(len(Joc.celuleGrid)):

                        if Joc.celuleGrid[np].collidepoint(pos):
                            linie = np//10
                            coloana = np % 10
                            ###############################
                            if(stare_curenta.tabla_joc.matr[np]) == Joc.GOL:   #punctul 6
                                #realizam mutarile una cate una
                                #se selecteaza prima mutare din cele doua
                                if jumate_piesa_plasata is False and stare_curenta.tabla_joc.mutare_valida(linie,coloana):
                                    stare_curenta.tabla_joc.matr[linie*10 +coloana] = Joc.JMIN
                                    linie_jumate_piesa = linie
                                    coloana_jumate_piesa = coloana
                                    jumate_piesa_plasata = True

                                elif jumate_piesa_plasata is True and stare_curenta.tabla_joc.pot_completa_piesa(linie_jumate_piesa,coloana_jumate_piesa,linie,coloana):
                                    stare_curenta.tabla_joc.matr[linie*10+coloana] = Joc.JMIN
                                    linie_jumate_piesa = -1
                                    coloana_jumate_piesa = -1
                                    jumate_piesa_plasata = False




                                # afisarea starii jocului in urma mutarii utilizatorului
                                print("\nTabla dupa mutarea jucatorului")
                                print(str(stare_curenta))

                                stare_curenta.tabla_joc.deseneaza_grid(coloana_marcaj=coloana)
                                # testez daca jocul a ajuns intr-o stare finala
                                # si afisez un mesaj corespunzator in caz ca da
                                if (afis_daca_final(stare_curenta)):
                                    break

                                # S-a realizat o mutare. Schimb jucatorul cu cel opus
                                stare_curenta.j_curent = Joc.jucator_opus(stare_curenta.j_curent)



        # --------------------------------
        else:  # jucatorul e JMAX (calculatorul)
            # Mutare calculator

            # preiau timpul in milisecunde de dinainte de mutare
            t_inainte = int(round(time.time() * 1000))
            if tip_algoritm == 'minimax':
                stare_actualizata = min_max(stare_curenta)
            else:  # tip_algoritm=="alphabeta"
                stare_actualizata = alpha_beta(-500, 500, stare_curenta)
            stare_curenta.tabla_joc = stare_actualizata.stare_aleasa.tabla_joc

            print("Tabla dupa mutarea calculatorului\n" + str(stare_curenta))

            # preiau timpul in milisecunde de dupa mutare
            t_dupa = int(round(time.time() * 1000))
            print("Calculatorul a \"gandit\" timp de " + str(t_dupa - t_inainte) + " milisecunde.")

            stare_curenta.tabla_joc.deseneaza_grid()
            if (afis_daca_final(stare_curenta)):
                print("Sfarsitul jocului")
                if(Joc.scor_x > Joc.scor_0):
                    print("A castigat X")
                else:
                    print("A castigat 0")

                break

            # S-a realizat o mutare. Schimb jucatorul cu cel opus
            stare_curenta.j_curent = Joc.jucator_opus(stare_curenta.j_curent)


if __name__ == "__main__":
    main()
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()