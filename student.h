/******************************************************************************
 * Projekt - Zaklady pocitacove grafiky - IZG
 * spanel@fit.vutbr.cz
 *
 * $Id:$
 */

#ifndef Student_H
#define Student_H

/******************************************************************************
 * Includes
 */

#include "render.h"

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
 * Studentsky renderer, do ktereho muzete dopisovat svuj kod
 */

/* Jadro vaseho rendereru */
typedef struct S_StudentRenderer
{
    /* Kopie default rendereru
     * Typicky trik jak implementovat "dedicnost" znamou z C++
     * Ukazatel na strukturu lze totiz pretypovat na ukazatel
     * na prvni prvek struktury a se strukturou S_StudentRenderer
     * tak pracovat jako s S_Renderer... */
    S_Renderer  base;

	S_RGBA *policka;

} S_StudentRenderer;


/******************************************************************************
 * Nasledujici fce budete nejspis muset re-implementovat podle sveho
 */

/* Funkce vytvori vas renderer a nainicializuje jej */
S_Renderer * studrenCreate();

/* Funkce korektne zrusi renderer a uvolni pamet */
void studrenRelease(S_Renderer **ppRenderer);

/* Nova fce pro rasterizaci trojuhelniku do frame bufferu
 * s podporou texturovani a interpolaci texturovacich souøadnic
 * Pozn.: neni nutné øešit perspektivní korekci textury
 * v1, v2, v3 - ukazatele na vrcholy trojuhelniku ve 3D pred projekci
 * n1, n2, n3 - ukazatele na normaly ve vrcholech ve 3D pred projekci
 * t1, t2, t3 - ukazatele na texturovaci souradnice vrcholu
 * x1, y1, ... - vrcholy trojuhelniku po projekci do roviny obrazovky */
void studrenDrawTriangle(S_Renderer *pRenderer,
                         S_Coords *v1, S_Coords *v2, S_Coords *v3,
                         S_Coords *n1, S_Coords *n2, S_Coords *n3,
                         S_Coords *t1, S_Coords *t2, S_Coords *t3,
                         int x1, int y1,
                         int x2, int y2,
                         int x3, int y3
                         );

/* Vykresli i-ty trojuhelnik modelu pomoci nove fce studrenDrawTriangle()
 * Pred vykreslenim aplikuje na vrcholy a normaly trojuhelniku
 * aktualne nastavene transformacni matice!
 * i - index trojuhelniku */
void studrenProjectTriangle(S_Renderer *pRenderer, S_Model *pModel, int i);

/* Vrací hodnotu v aktuálnì nastavené textuøe na zadaných
 * texturovacích souøadnicích u, v
 * Pro urèení hodnoty používá bilineární interpolaci
 * u, v - texturovací souøadnice v intervalu 0..1, který odpovídá šíøce/výšce textury */
S_RGBA studrenTextureValue(S_StudentRenderer * pRenderer, double u, double v);


/******************************************************************************
 * Deklarace fci, ktere jsou umistene ve student.c a pouzivaji se v main.cc
 */

/* Callback funkce volana pri startu aplikace
 * Doplnte automaticke vygenerovani a prirazeni texturovacich souradnic
 * vrcholum modelu s vyuzitim mapovani na kouli
 * Muzete predpokladat, ze model je umisten v pocatku souradneho systemu
 * a posunuti neni treba resit */
void onInit(S_Renderer *pRenderer, S_Model *pModel);


#ifdef __cplusplus
}
#endif

#endif /* Student_H */

/*****************************************************************************/
/*****************************************************************************/
