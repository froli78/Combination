#ifndef COMBI_LIB_H_INCLUDED
#define COMBI_LIB_H_INCLUDED

 /** Combinations - (2019.12.20)
 *
 * combinations.cpp
 * Copyright (C) 2019 - Filip Roland - git@filiproland.hu
 *
 * Combinations is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, version 3 of the License.
 *
 * Combinations is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define _UNICODE
#define UNICODE

#include <string>
#include <fstream> //Filekezelés - írás és olvasás vegyesen
#include <future> //Aszinkron szál futtatás

#include <limits>

#include <binbigint.h>

typedef unsigned long long int ulli;


namespace combinatorics
{

    const ulli ulliLimit = (ulli)std::numeric_limits<ulli>::max;


    //Faktoriális számítás (mivel az ulli = 2^64 = 1.8E+19 ezért ezzel max 20!-ig lehet számolni (mivel 20! = 2,4E+18 de már 21! = 5,1E+19...)
    ulli factorial
    (
        ulli f
    )
    {
        return (f == 1 || f == 0) ? 1 : factorial(f - 1) * f;
    }

    //Bigint faktoriális függvény
    BinBigInt biFactorial
    (
        const BinBigInt & f    //A faktoriális száma
    )
    {
        BinBigInt r;   //Az eredményt tartalmazni hivatott BinBigInt
        BinBigInt::biFactorial(f.v, f.sigNeg, r.v, r.sigNeg);
        return r;
    }


    //Bigint produktum függvény
    BinBigInt biPosProductOfSequence
    (
        const BinBigInt & i,     //A produktum kezdő száma
        const BinBigInt & n      //A produktum befejező száma
    )
    {
        BinBigInt r;   //Az eredményt tartalmazni hivatott BinBigInt
        BinBigInt::biPosProductOfSequence(i.v, i.sigNeg, n.v, n.sigNeg, r.v, r.sigNeg);
        return r;
    }

    //Bigint produktum függvény
    BinBigInt biProductOfSequence
    (
        BinBigInt i,     //A produktum kezdő száma
        BinBigInt n      //A produktum befejező száma
    )
    {
        BinBigInt r;   //Az eredményt tartalmazni hivatott BinBigInt
        if (i.sigNeg == 1 || n.sigNeg == 1)
        {
            if (i.sigNeg == 1 && n.sigNeg == 1) //ha i és n is negatív
            {
                i = i.biAbs();
                n = n.biAbs();
                //i.biAbs(i);
                //n.biAbs(n);

                BinBigInt::biPosProductOfSequence(i.v, i.sigNeg, n.v, n.sigNeg, r.v, r.sigNeg);

                //Az elemszám páros vagy páratlan voltának megállapítása
                BinBigInt one((ulli)1);
                BinBigInt d = n-i+one;
                d = d.biAbs();
                //d.biAbs(d);
                BinBigInt two((ulli)2);
                BinBigInt odd = d % two;
                if (!BinBigInt::biIsZero(odd)) //ha az elemszám páratlan
                {
                    r.BinBigInt::biNegative(); //Az eredmény negatívra váltása
                }
            }
            else //Ha csak az egyik negatív
            {
                if (i <= n) //Ha i <= n
                {
                    i = i.biAbs();
                    //i.biAbs(i);
                    BinBigInt one((ulli)1);

                    BinBigInt rNegPart;

                    BinBigInt::biPosProductOfSequence(one.v, one.sigNeg, i.v, i.sigNeg, rNegPart.v, rNegPart.sigNeg);

                    //Az elemszám páros vagy páratlan voltának megállapítása (itt az elemszám = i-vel)
                    BinBigInt two((ulli)2);
                    BinBigInt odd = i % two;
                    if (!BinBigInt::biIsZero(odd)) //ha az elemszám páratlan
                    {
                        rNegPart.BinBigInt::biNegative(); //Az eredmény negatívra váltása
                    }

                    BinBigInt::biPosProductOfSequence(one.v, one.sigNeg, n.v, n.sigNeg, r.v, r.sigNeg);

                    r = r * rNegPart;
                }
                else
                {
                    n = n.biAbs();
                    //n.biAbs(n);
                    BinBigInt one((ulli)1);

                    BinBigInt rNegPart;

                    BinBigInt::biPosProductOfSequence(one.v, one.sigNeg, n.v, n.sigNeg, rNegPart.v, rNegPart.sigNeg);

                    //Az elemszám páros vagy páratlan voltának megállapítása (itt az elemszám = i-vel)
                    BinBigInt two((ulli)2);
                    BinBigInt odd = n % two;
                    if (!BinBigInt::biIsZero(odd)) //ha az elemszám páratlan
                    {
                        rNegPart.sigNeg = 1; //Az eredmény előjele negatív
                    }

                    BinBigInt::biPosProductOfSequence(one.v, one.sigNeg, i.v, i.sigNeg, r.v, r.sigNeg);

                    r = r * rNegPart;
                }
            }
        }
        else //Ha i és n is pozitív
        {
                    BinBigInt::biPosProductOfSequence(i.v, i.sigNeg, n.v, n.sigNeg, r.v, r.sigNeg);
        }

        return r;
    }

    //Bigint nCk (n alatt a k) függvény
    BinBigInt binCk
    (
        const BinBigInt & n,     //n
        const BinBigInt & k      //k
    )
    {
        BinBigInt one((ulli)1);
        BinBigInt nmkpone = n-k+one; //BinBigInt (n-k+1)

        BinBigInt nprod;
        BinBigInt::biPosProductOfSequence(nmkpone, n, nprod); //n = ∏ (n-k+1) to n = (n-k+1)*(n-k+2)*..*(n-1)*n = n!/(n-k)!

        BinBigInt kfact;
        BinBigInt::biFactorial(k, kfact); //k!

        return (nprod/kfact); //nCk = n!/(n-k)!/k!
    }

// TODO (FR#1#): Kellene egy fájlból betöltő eljárás is...
// TODO (FR#1#): Szükséges lenne még egy az 1toN alap kombinációkon alapuló sima elem egyetzetést végző és lecserélő függvény
//illetve egy tetszőleges elemeket tartalmazó tömbbel működőközvetlenül egyedi elem kombinációkat előállító verzió is...

    template <class T>
    class NonrepCombination
    {

            //Tagváltozók
            ulli n;
            ulli k;
            T * headerArray;                //T típusú pointer dekarálás dinamikus memória foglaláshoz fejléc számára
            T * combinationDataArray;       //T típusú pointer dekarálás dinamikus memória foglaláshoz a kombinációs rekordok számára
            unsigned char headerSize = 25;  //Fejléc méret (fájlba íráshoz szüskéges)
            ulli dataSize = 0;              //A kombinációs rekord tömb méretének tárolására szolgáló tag változó (fájlba íráshoz szüskéges)
            bool generated = false;         //A kombináció létrehozásnak megtörténtét jelző flag
            bool saved = false;             //A kombináció fájlba mentésének megtörténtét jelző flag

        public:

            //KONSTRUSKTOROK, DESTRUKTOR

            NonrepCombination(ulli n = 0, ulli k = 0) //:n=n, k=k {}
            {
                this->n=n;
                this->k=k;
            }

            ~NonrepCombination()
            {
                delete [] combinationDataArray; //Dinamikus memória felszabadítása
                delete [] headerArray; //Dinamikus memória felszabadítása
            }

        private:

            //Tagfüggvények

            //Az Ismétlésnélküli kombinációkat elõállító függvény által használt rekurzív függvény ami visszaadja a hátulról értelmezett soron következő lépésben megnövelendő kombinációs rekordelem pozícióját
            //template <class T>
            //void NonrepCombination<T>::combLimitChkFromEnd
            void combLimitChkFromEnd
            (
                ulli n, //n
                ulli k, //k
                ulli &p, //az aktuálisan növelendő hátulról értelmezett inverz pozíció index kiinduló értéke (induláskor = 0 kell legyen!) és a végső érték visszadására szolgáló referencia
                T * crp //Az aktuálisan vizsgált k elemű kombinációs rekord pointere
            )
            {
                if (ulli(crp[k-1-p]) == n-p && p < k) //Ha az aktuális p szerinti pozícióban található érték már elérte az arra a pozícióra értelmezett maximum értéket (=n-p)
                {
                    p=p+1; //A p növelése eggyel (egy pozícióval elõrébb [a k elemű kombinációs rekord eleje felé] lép)
                    NonrepCombination<T>::combLimitChkFromEnd(n, k, p, crp); //Az elõrébb lévõ pozíció rekurzív ellenõrzése (Ha valahol még nem teljesül a feltétel kilép)
                    //NonrepCombination<T>::combLimitChkFromEnd<T>(n, k, p, crp); //Az elõrébb lévõ pozíció rekurzív ellenõrzése (Ha valahol még nem teljesül a feltétel kilép)
                }
            }

            //template <class T>
            //void NonrepCombination<T>::combLimitChkFromBegin
            void combLimitChkFromBegin
            (
                ulli n, //n
                ulli k, //k
                ulli &p, //az aktuálisan növelendő hátulról értelmezett inverz pozíció index kiinduló értéke (induláskor = 0 kell legyen!) és a végső érték visszadására szolgáló referencia
                T * crp //Az aktuálisan vizsgált k elemű kombinációs rekord pointere
            )
            {
                if (crp[p] == crp[p-1]+1 && p > 0) //Ha az aktuális p szerinti pozícióban található érték egyel nagyobb mint a tőle balra álló érték és p > 0
                {
                    p=p-1; //A p csökkentése eggyel (egy pozícióval elõrébb [a k elemű kombinációs rekord eleje felé] lép)
                    combLimitChkFromBegin(n, k, p, crp); //Az elõrébb lévõ pozíció rekurzív ellenõrzése (Ha valahol még nem teljesül a feltétel kilép)
                    //combLimitChkFromBegin<T>(n, k, p, crp); //Az elõrébb lévõ pozíció rekurzív ellenõrzése (Ha valahol még nem teljesül a feltétel kilép)
                }
            }

            //A kombinációk előállítása és memóriába helyezése (Előrefelé verzió)
            //template <class T>
            //ulli NonrepCombination<T>::forwardCombGeneration
            ulli forwardCombGeneration
            (
                ulli n, //n
                ulli k, //k
                T * combinationDataArray,  //A kombinációs rekordokat tároló memória tömb
                ulli forwardRecordNr = ulliLimit //A generálandó rekordok száma ha nincs megadva akkor ulli_MAX
            )
            {
                ulli createdRecordNr;

                T * crp;

                try
                {
                    crp = new T [k]{0}; //Egy kombinációs rekord tárolására alkalmas tömb az egyes kombinációk előállítására (Minden elem 0-ra inicializálva.)
                }
                catch (const std::bad_alloc e)
                {
                        printf("Exception: ");
                        printf(e.what());
                        printf("\n");
                    //MessageBoxA(NULL, e.what(), "Exception", MB_OK);

                        printf("AsyncThread: Nincs elég memória. A szükséges memória mérete = %llu byte!\n", (k)*sizeof(T));
                    /*
                    std::wstring stro1, stro2;
                    stro1 = L"Nincs elég memória. A szükséges memória mérete = ";
                    stro2 = std::to_wstring((k)*sizeof(T));
                    stro1.append(stro2);
                    stro2 = L" byte!";
                    stro1.append(stro2);
                    MessageBox(NULL, stro1.c_str(), L"Exception", MB_OK);
                    */
                    //MessageBox(NULL, L"Nincs elég memória", L"Exception", MB_OK);
        // TODO (FR#1#): Ez esetben vagy közvetlen háttértárra írással történő megoldás alkalmazása  vagy a rendelkezésre álló memóriának megfelelő részenkénti előállítás (lapozás) és részenként kiírás lehetséges...

                    return 0;
                }
                 catch (...)
                {
                    printf("Exception: Általános hiba\n");
                    //MessageBox(NULL, L"Általános hiba", L"Exception", MB_OK);
                    return 0;
                }

                //Sikeres memória foglalás visszajelzése (csak debug infó jelleggel)
                //printf("AsyncThread: A memória foglalás sikeres volt. A lefoglalt memória mérete =  %llu byte!\n", (k)*sizeof(T));

                //a kombinációs rekord inicializálása {1, 2, 3, ... , k-2, k-1, k}
                for (ulli i=0; i < k; i++)
                {
                    crp[i] = i+1; //a kombinációs rekord inicializálása {1, 2, 3, ... , k-2, k-1, k}
                }

                ulli p=0; //Az aktuális k elemű kombinációs rekordon belül hátulról az aktuálisan megnövelendő pozíció inverz indexe hátulról 0-tól számozva kezdődően
                createdRecordNr = 0; //Tényleges rekord számlálás (1-től kezdődően értelmezve) és egyben a rekord pozíció indexálás (0-tól kezdődően értelmezve) változója
                while (p < k && createdRecordNr < forwardRecordNr) //Amíg az aktuális rekordon beüli megnövelendő inverz pozíció index legfeljebb a legelső (k-1 inverz indexű) pozícióra mutat
                {

                    //A soron következő azon kombinációs rekordok előállítása amelyekben csak az utolsó elem változik (ezek egyszerű for cikussal előállíthatók)
                    while (crp[k-1] <= n && createdRecordNr < forwardRecordNr) //Amíg az aktuális kombináció utolsó eleme nem nagyobb mint n
                    {
                        //Az aktuális kombinációs rekord memóriába helyezése
                        for (ulli c = 0; c < k; c++) //k× (0-->k-1 indexig)
                        {
                            combinationDataArray[createdRecordNr*k+c] = crp[c]; //Az aktuális rekord elem a memóriában az aktuális rekordnak megfelelő pozícióba, ahol a (createdRecordNr*k) az aktuális rekordnak megfelelő pozíció a memória tömbben
                        }
                        createdRecordNr++; //Tényleges rekord számlálás
                        crp[k-1]++; //Az utolsó rekordelem növelése egyel (a soron következő olyan rekord amiben csak az utolsó elem változik)
                    }
                    crp[k-1]--; //Az utolsó rekod elem vissza léptetése n-re
                    p = 0; //p a hátulról az aktuálisan megnövelendő pozíció inverz indexének minden cikusban történő nullázása (a k elemű kombináció utolsó elemére visszaállítása, fontos, hogy a végéről kezdje ellenőrizni!)
                    combLimitChkFromEnd(n, k, p, crp); //A hátulról az aktuálisan megnövelendő pozíció inverz indexe meghatározása
                    crp[k-1-p]++; //A hátulról az aktuálisan soron következő megnövelendő pozícióban lévő érték megnövelése
                    for (ulli r = 1; r<=p; r++) //A k elemű kombinációs rekordban az aktuális pozíciótól az utolsó pozícióig 1-esével növekvő elemek beállítása
                    {
                        crp[k-1-p+r] = crp[k-1-p+r-1]+1; //Az aktuális pozíciótól a kombinációs rekord tömben sorban következő értékek mindegyikének az előtte lévőnél eggyel nagyobbra álíltása (a következő cikus kiinduló értékének beálíltása)
                    }
                }

                delete [] crp; //Dinamikus memória felszabadítása

                //printf("Debug: Előröl verzió: kombináció kész\n");
                //MessageBox(NULL, L"Kombináció kész", L"Debug", MB_OK);

                return createdRecordNr;
            }

            //A kombinációk előállítása és memóriába helyezése (Hátrafelé verzió)
            //template <class T>
            //ulli NonrepCombination<T>::backwardCombGeneration
            ulli backwardCombGeneration
            (
                ulli n, //n
                ulli k, //k
                T * combinationDataArray,  //A kombinációs rekordokat tároló memória tömb
                ulli recordNr, //A generálandó rekordok teljes száma (tehát a forward elemszámmal együttes elemszám, a tömb végének meghatározásához szüskéges)
                ulli backwardRecordNr = ulliLimit //A generálandó rekordok száma ha nincs megadva akkor ulli_MAX
            )
            {
                ulli createdRecordNr;

                T * crp;

                try
                {
                    crp = new T [k]{0}; //Egy kombinációs rekord tárolására alkalmas tömb az egyes kombinációk előállítására (Minden elem 0-ra inicializálva.)
                }
                catch (const std::bad_alloc e)
                {
                        printf("Exception: ");
                        printf(e.what());
                        printf("\n");
                    //MessageBoxA(NULL, e.what(), "Exception", MB_OK);
                        printf("AsyncThread: Nincs elég memória. A szükséges memória mérete = %llu byte!\n", (k)*sizeof(T));
                    /*
                    std::wstring stro1, stro2;
                    stro1 = L"Nincs elég memória. A szükséges memória mérete = ";
                    stro2 = std::to_wstring((k)*sizeof(T));
                    stro1.append(stro2);
                    stro2 = L" byte!";
                    stro1.append(stro2);
                    MessageBox(NULL, stro1.c_str(), L"Exception", MB_OK);
                    */
                    //MessageBox(NULL, L"Nincs elég memória", L"Exception", MB_OK);
        // TODO (FR#1#): Ez esetben vagy közvetlen háttértárra írással történő megoldás alkalmazása  vagy a rendelkezésre álló memóriának megfelelő részenkénti előállítás (lapozás) és részenként kiírás lehetséges...

                    return 0;
                }
                 catch (...)
                {
                    printf("Exception: Általános hiba\n");
                    //MessageBox(NULL, L"Általános hiba", L"Exception", MB_OK);
                    return 0;
                }

                //Sikeres memória foglalás visszajelzése (csak debug infó jelleggel)
                //printf("AsyncThread: A memória foglalás sikeres volt. A lefoglalt memória mérete =  %llu byte!\n", (k)*sizeof(T));

                //a kombinációs rekord inicializálása {n-k+1, n-k+2, n-k+3, ..., n-2, n-1 ,n} (n-től lefelé k db elem)
                for (ulli i=0; i < k; i++)
                {
                    crp[i] = i+n-k+1; //a kombinációs rekord inicializálása {n-k+1, n-k+2, n-k+3, ..., n-2, n-1 ,n} (n-től lefelé k db elem)
                }

                ulli p = k-1; //Az aktuális k elemű kombinációs rekordon belül elölről az aktuálisan csökkentendő pozíció indexe 0-tól számozva
                createdRecordNr = 0; //Tényleges rekord számlálás (1-től kezdődően értelmezve) és egyben a rekord pozíció indexálás (0-tól kezdődően értelmezve) változója
                while (crp[p] >= p+1 && createdRecordNr < backwardRecordNr) //Amíg az aktuális rekordon beüli csökkentendő pozíció index legfeljebb az utolsó (k-1 indexű) pozícióra mutat
                {
                    //A soron következő azon kombinációs rekordok előállítása amelyekben csak az utolsó elem változik (ezek egyszerű for cikussal előállíthatók)
                    while (crp[k-1] > crp[k-2]  && createdRecordNr < backwardRecordNr) //Amíg az aktuális kombináció utolsó eleme nagyobb mint az utolsó előtti elem
                    {
                        //Az aktuális kombinációs rekord memóriába helyezése (hátulról előre, de a kombinációs rekord elemeinek normál sorendjében)
                        for (ulli c = 0; c < k; c++) //k× (0-->k-1 indexig)
                        {
                            combinationDataArray[(recordNr-1-createdRecordNr)*k+c] = crp[c]; //Az aktuális rekord elem a memóriában az aktuális rekordnak megfelelő pozícióba, ahol a ((recordNr-1-createdRecordNr)*k) az aktuális rekordnak megfelelő pozíció (hátulról) a memória tömbben
                        }
                        createdRecordNr++; //Tényleges rekord számlálás
                        crp[k-1]--; //Az utolsó rekordelem csökkentése egyel (a soron következő olyan rekord amiben csak az utolsó elem változik)
                    }
                    crp[k-1]++; //Az utolsó rekod elem vissza léptetése egyel

                    p = k-1; //p az elölről az aktuálisan csökkentendő pozíció inverz indexének minden cikusban történő k-1-re állítása (a k elemű kombináció utolsó elemére visszaállítása, fontos, hogy a végéről kezdje ellenőrizni!)
                    combLimitChkFromBegin(n, k, p, crp); //A hátulról az aktuálisan csökkentendő pozíció indexe meghatározása
                    crp[p]--; //Az aktuálisan soron következő csökkentendő pozícióban lévő érték csökkentése
                    for (ulli r = 0; r<k-p-1; r++) //A k elemű kombinációs rekordban az aktuális pozíciótól az utolsó pozícióig hátulról n-től 1-esével csökkenő elemek beállítása
                    {
                        crp[k-1-r] = n-r; //A k elemű kombinációs rekordban az aktuális pozíciótól az utolsó pozícióig hátulról n-től 1-esével csökkenő elemek beállítása
                    }
                }

                delete [] crp; //Dinamikus memória felszabadítása

                //printf("Debug: Hátulról verzió: kombináció kész\n");
                //MessageBox(NULL, L"Kombináció kész", L"Debug", MB_OK);

                return createdRecordNr;
            }

        public:

            //Ismétlésnélküli kombinációkat előállító függvény
            //template <class T>
            //ulli NonrepCombination<T>::generate1ToNNonrepCombination
            ulli generate1ToNNonrepCombination()
            {
                if (n == 0 || k == 0) throw std::invalid_argument("Exception: None of the parameters can be zero!");

                //printf("nonrepCombination() függvény elindult\n");
                //MessageBox(NULL, L"MessageBox függvényből.", L"Próbass", MB_OK);

                //Az ismétlés nélküli kombináció elemszámának számítása
                //Iinicializálás
                ulli ulliMaxV, kf, nmktnf, i, r, recordNr;
                recordNr = 0;

                //Kombináció szám számítása

                std::string str;

                BinBigInt bin((ulli)n); //BinBigInt n
                BinBigInt bik((ulli)k); //BinBigInt k
                BinBigInt birecordNr; //A kombinációk számának tárolására létrehozott BinBigInt változó

                //birecordNr = combinatorics::binCk(n, k);
                birecordNr = binCk(bin, bik);

                //Annak ellenőrzése, hogy a kombináció rekordszám számítás eredeménye belefér-e egy unsigned lon long int változóba...
                //Ha nem jeleneg kivételt dob de ez nem feltétlenül szükséges ha erre a változóra a további kódban is ezt a BinBigInt változót használnám csak ez még kód átalakítást igényelne
                BinBigInt biulliMax(BinBigInt::vUlliLimit, 0);
                if (birecordNr > biulliMax)
                {
                    //printf("A kombináció rekordszáma meghaladja az ulli változó limitjét.\n");
                    throw std::range_error("A kombináció rekordszáma meghaladja az ulli változó limitjét.");
                }

                //Konverzió ullire
                recordNr = (ulli)*birecordNr.data();

                //printf("recordNr ulli konverzió utáni ereménye = %llu\n", recordNr);

                //Annak ellenőrzése, hogy a kombináció elemszáma belefér-e egy unsigned lon long int változóba...
                //Ha nem jeleneg kivételt dob de ez nem feltétlenül szükséges ha a további kódban a feladatot részekre bontva végezném el...
                if (birecordNr*k > BinBigInt::ulliLimit)
                {
                    //printf("A kombináció elemszáma meghaladja az ulli változó limitjét.\n");
                    throw std::range_error("A kombináció elemszáma meghaladja az ulli változó limitjét.");
                }

                //KOMBINÁCIÓK ELŐÁLLÍTÁSA MEMÓRIÁBA

                //Memória foglalás
                try
                {
                    headerArray = new T [headerSize]{0}; //A kombinációs fájl header adatait tartalmazó 25 byteos tömb
                    combinationDataArray = new T [recordNr*k]{0}; //Az eredmény tárolására szolgáló memória foglalása

                }
                catch (const std::bad_alloc e)
                {
                        printf("Exception: ");
                        printf(e.what());
                        printf("\n");
                    //MessageBoxA(NULL, e.what(), "Exception", MB_OK);
                        printf("Nincs elég memória. A szükséges memória mérete = %llu byte!\n", (recordNr*k+headerSize)*sizeof(T));
                    /*
                    std::wstring stro1, stro2;
                    stro1 = L"Nincs elég memória. A szükséges memória mérete = ";
                    stro2 = std::to_wstring((recordNr*k+headerSize)*sizeof(T));
                    stro1.append(stro2);
                    stro2 = L" byte!";
                    stro1.append(stro2);
                    MessageBox(NULL, stro1.c_str(), L"Exception", MB_OK);
                    */
                    //MessageBox(NULL, L"Nincs elég memória", L"Exception", MB_OK);
        // TODO (FR#1#): Ez esetben vagy közvetlen háttértárra írással történő megoldás alkalmazása  vagy a rendelkezésre álló memóriának megfelelő részenkénti előállítás (lapozás) és részenként kiírás lehetséges...

                    return 0;
                }
                catch (...)
                {
                    printf("Exception: Általános hiba\n");
                    //MessageBox(NULL, L"Általános hiba", L"Exception", MB_OK);
                    return 0;
                }

                //Sikeres memória foglalás visszajelzése (csak debug infó jelleggel)
                //printf("A memória foglalás sikeres volt. A lefoglalt memória mérete =  %llu byte!\n", (recordNr*k+headerSize)*sizeof(T));
                /*
                std::wstring stro1, stro2;
                stro1 = L"A memória foglalás sikeres volt. A lefoglalt memória mérete = ";
                stro2 = std::to_wstring((recordNr*k)*sizeof(T));
                stro1.append(stro2);
                stro2 = L" byte!";
                stro1.append(stro2);
                MessageBox(NULL, stro1.c_str(), L"Memória", MB_OK);
                */

                ulli forwardRecordNr = recordNr/2+recordNr%2; //Az előrefelé generálandó kombinációk száma
                printf("forwardRecordNr =  %llu\n", forwardRecordNr);

                ulli backwardRecordNr = recordNr/2;  //A hárafelé generálandó kombinációk kombinációk
                printf("backwardRecordNr =  %llu\n", backwardRecordNr);


                //Aszinkron szál indítása
                std::future<ulli> fCreatedForwardRecordNr = std::async(&NonrepCombination<T>::forwardCombGeneration, this, n, k, combinationDataArray, forwardRecordNr);

                //Main thread
                ulli createdBackwardRecordNr = backwardCombGeneration(n, k, combinationDataArray, recordNr,  backwardRecordNr);

                ulli createdForwardRecordNr = 0;

                createdForwardRecordNr = fCreatedForwardRecordNr.get(); //Aszinkron szál érték visszaadása (bevárása)


                printf("createdForwardRecordNr = %llu\n", createdForwardRecordNr);
                printf("createdBackwardRecordNr = %llu\n", createdBackwardRecordNr);

                ulli createdRecordNr = 0;
                createdRecordNr = createdForwardRecordNr + createdBackwardRecordNr;
                printf("createdRecordNr = %llu\n", createdRecordNr);

                //A tényleges eredmény kombinációk és a számítható kombinációk számának összehasonítása ellenőrzés céljából
                if (createdRecordNr != recordNr)
                {
                    createdRecordNr = 0;
                    printf("Hiba: A tényleges kombinációk száma nem egyezik a számított értékkel!\n");
                    //MessageBox(NULL, L"A tényleges kombinációk száma nem egyezik a számított értékkel!", L"Hiba:", MB_OK);
                }


                //Header előállítása
                //KELL 25 byte az elejére: 1-8. n (qword), 9-16. k (qword), 17-24. nCk = rekodszám (qword), 25. elem méret [1-byte, 2-word, 4-dword, 8-qword] (byte)
                ulli temp = n;
                for (i=0; i <= 7; i++)
                {
                    headerArray[i] = temp; //Az egy kombinációban (kombinációs rekordban) található elemek száma = k
                    temp = temp >> 8; //SHIFT right 8 (következő byte az alsó bytba)
                }
                temp = k;
                for (i=8; i <= 15; i++)
                {
                    headerArray[i] = temp; //Az egy kombinációban (kombinációs rekordban) található elemek száma = k
                    temp = temp >> 8; //SHIFT right 8 (következő byte az alsó bytba)
                }
                temp = recordNr;
                for (i=16; i <= 23; i++)
                {
                    headerArray[i] = temp; //A kombinációk (kombinációs rekordok) száma = n alatt k
                    temp = temp >> 8; //SHIFT right 8 (következő byte az alsó bytba)
                }
                headerArray[24] = char(sizeof(T)); //Egy adat érték byteokban számított értéke (1-byte [char], 2-word [short int / int], 4-dword [int / long int], 8-qword [long int / long long int])

                dataSize = recordNr*k*sizeof(T); //A kombinációk adatmérete

                generated = true; //A kombináció létrehozásnak megtörténtét jelző flag beállítása

                return createdRecordNr;
            }

            std::string writeToFile
            (
                std::string fileDirectoryPath
            )
            {
                if (generated)
                {
                    //Memória alapú eremény Fájlba írása
                    //Kiementi fájl teljes elérési útjának és fájlnevének létrehozása
                    std::string stro1, stro2;
                    stro1 = fileDirectoryPath;
                    stro2 = std::to_string(n);
                    stro1.append(stro2);
                    stro2 = "c";
                    stro1.append(stro2);
                    stro2 = std::to_string(k);
                    stro1.append(stro2);
                    stro2 = ".cmb";
                    stro1.append(stro2);
                    fileDirectoryPath = stro1;

                    //printf("A kimeneti fájl teljes elérési útja: %s\n", fileDirectoryPath.c_str());

                /*
                    //Kiementi fájl teljes elérési útjának és fájlnevének létrehozása
                    stro1 = fileDirectoryPath;
                    stro2 = std::to_wstring(n);
                    stro1.append(stro2);
                    stro2 = L"a";
                    stro1.append(stro2);
                    stro2 = std::to_wstring(k);
                    stro1.append(stro2);
                    stro2 = L".cmb";
                    stro1.append(stro2);
                    fileDirectoryPath = stro1;
                */
                /*
                    stro1 = L"A kimeneti fájl teljes elérési útja:";
                    stro2 = fileDirectoryPath;
                    stro1.append(stro2);
                    MessageBox(NULL, stro1.c_str(), L"Infó:", MB_OK);
                 //*/

                    //Kimeneti fájl létrehozása
                    std::fstream f; //fájl változó dekarálása
                    f.open(fileDirectoryPath, std::ios::out | std::ios::binary);
                    f.write(reinterpret_cast<const char *>(&headerArray[0]), headerSize);
                    f.write(reinterpret_cast<const char *>(&combinationDataArray[0]), dataSize);

                    f.close();
                    //printf("Debug: Fájlbaírás után\n");
                    //MessageBox(NULL, L"Fájlbaírás után", L"Debug", MB_OK);
                }
                else if (saved)
                {
                    printf("Exception: The combination saved yet. Nothing to do.");
                }
                else
                {
                    throw "Exception: The combination not generated yet. Nothing to write to file.";
                }

                return fileDirectoryPath;
            }
    };
}

#endif // COMBI_LIB_H_INCLUDED
