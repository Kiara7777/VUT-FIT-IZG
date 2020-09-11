/******************************************************************************
 * Projekt - Zaklady pocitacove grafiky - IZG
 * spanel@fit.vutbr.cz
 *
 * $Id: student.c 337 2014-02-25 06:52:49Z spanel $
 */

#include "student.h"
#include "transform.h"
#include "bmp.h"

#include <memory.h>
#include <math.h>
#include <stdbool.h>


/******************************************************************************
 * Globalni promenne a konstanty
 */

/* rozmer textury */
const int       TEXTURE_SIZE    = 512;

/* pocet policek sachovnice */
const int       NUM_OF_TILES    = 16;

/* barva poli */
const S_RGBA    BLACK_TILE      = { 75, 75, 75 };
const S_RGBA    WHITE_TILE      = { 255, 255, 255 };


/*****************************************************************************
 * Funkce vytvori vas renderer a nainicializuje jej
 */

S_Renderer * studrenCreate()
{
    S_StudentRenderer * renderer = (S_StudentRenderer *)malloc(sizeof(S_StudentRenderer));
    IZG_CHECK(renderer, "Cannot allocate enough memory");

    /* inicializace default rendereru */
    renInit(&renderer->base);

    /* nastaveni ukazatelu na upravene funkce */
	renderer->base.releaseFunc = studrenRelease;
	renderer->base.projectTriangleFunc = studrenProjectTriangle;

    /* inicializace nove pridanych casti */
    /* ??? */
	renderer->policka = malloc(TEXTURE_SIZE * TEXTURE_SIZE * sizeof(S_RGBA));
	IZG_CHECK(renderer, "Cannot allocate enough memory");

	int pixel_sloupec = 0;
	int pixel_radek = 0;
	bool white = true;


	//vytvoreni textury
	S_RGBA color = WHITE_TILE;

	for (int i = 0; i < TEXTURE_SIZE; i++)
	{
		for (int j = 0; j < TEXTURE_SIZE; j++)
		{
			if (pixel_sloupec < TEXTURE_SIZE / NUM_OF_TILES)
			{
				pixel_sloupec++;
			}
			else
			{
				if (pixel_radek >= TEXTURE_SIZE / NUM_OF_TILES)
				{
					pixel_sloupec = 1;
					pixel_radek = 0;
					
				}
				else
				{
					if (white)
					{
						white = false;
						color = BLACK_TILE;
						pixel_sloupec = 1;
					}
					else
					{
						white = true;
						color = WHITE_TILE;
						pixel_sloupec = 1;
					}
				}
			}
			renderer->policka[i * TEXTURE_SIZE + j] = color;
		
		}
		pixel_radek++;
	}
	
	//saveBitmap("textura", renderer->policka, TEXTURE_SIZE, TEXTURE_SIZE);
	return (S_Renderer *)renderer;
}

/*****************************************************************************
 * Funkce korektne zrusi renderer a uvolni pamet
 */

void studrenRelease(S_Renderer **ppRenderer)
{
    S_StudentRenderer * renderer;

    if( ppRenderer && *ppRenderer )
    {
        renderer = (S_StudentRenderer *)(*ppRenderer);

        /* pripadne uvolneni pameti */
         free(renderer->policka);
        /* fce default rendereru */
        renRelease(ppRenderer);
    }
}

/******************************************************************************
 * Nova fce pro rasterizaci trojuhelniku do frame bufferu
 * s podporou texturovani a interpolaci texturovacich souøadnic
 * Pozn.: neni nutné øe?it perspektivní korekci textury
 * v1, v2, v3 - ukazatele na vrcholy trojuhelniku ve 3D pred projekci
 * n1, n2, n3 - ukazatele na normaly ve vrcholech ve 3D pred projekci
 * t1, t2, t3 - ukazatele na texturovaci souradnice vrcholu
 * x1, y1, ... - vrcholy trojuhelniku po projekci do roviny obrazovky
 */

void studrenDrawTriangle(S_Renderer *pRenderer,
                         S_Coords *v1, S_Coords *v2, S_Coords *v3,
                         S_Coords *n1, S_Coords *n2, S_Coords *n3,
                         S_Coords *t1, S_Coords *t2, S_Coords *t3,
                         int x1, int y1,
                         int x2, int y2,
                         int x3, int y3
                         )
{
	int         minx, miny, maxx, maxy;
	int         a1, a2, a3, b1, b2, b3, c1, c2, c3;
	int         s1, s2, s3;
	int         x, y, e1, e2, e3;
	double      alpha, beta, w1, w2, w3, z, u, v;
	S_RGBA      col1, col2, col3, color;
	S_RGBA		textura;
	double		R, G, B;
	S_StudentRenderer *renderer = (S_StudentRenderer *)pRenderer;

	IZG_ASSERT(pRenderer && v1 && v2 && v3 && n1 && n2 && n3);

	/* vypocet barev ve vrcholech */
	col1 = pRenderer->calcReflectanceFunc(pRenderer, v1, n1);
	col2 = pRenderer->calcReflectanceFunc(pRenderer, v2, n2);
	col3 = pRenderer->calcReflectanceFunc(pRenderer, v3, n3);

	/* obalka trojuhleniku */
	minx = MIN(x1, MIN(x2, x3));
	maxx = MAX(x1, MAX(x2, x3));
	miny = MIN(y1, MIN(y2, y3));
	maxy = MAX(y1, MAX(y2, y3));

	/* oriznuti podle rozmeru okna */
	miny = MAX(miny, 0);
	maxy = MIN(maxy, pRenderer->frame_h - 1);
	minx = MAX(minx, 0);
	maxx = MIN(maxx, pRenderer->frame_w - 1);

	/* Pineduv alg. rasterizace troj.
	hranova fce je obecna rovnice primky Ax + By + C = 0
	primku prochazejici body (x1, y1) a (x2, y2) urcime jako
	(y1 - y2)x + (x2 - x1)y + x1y2 - x2y1 = 0 */

	/* normala primek - vektor kolmy k vektoru mezi dvema vrcholy, tedy (-dy, dx) */
	a1 = y1 - y2;
	a2 = y2 - y3;
	a3 = y3 - y1;
	b1 = x2 - x1;
	b2 = x3 - x2;
	b3 = x1 - x3;

	/* koeficient C */
	c1 = x1 * y2 - x2 * y1;
	c2 = x2 * y3 - x3 * y2;
	c3 = x3 * y1 - x1 * y3;

	/* vypocet hranove fce (vzdalenost od primky) pro protejsi body */
	s1 = a1 * x3 + b1 * y3 + c1;
	s2 = a2 * x1 + b2 * y1 + c2;
	s3 = a3 * x2 + b3 * y2 + c3;

	/* normalizace, aby vzdalenost od primky byla kladna uvnitr trojuhelniku */
	if (s1 < 0)
	{
		a1 *= -1;
		b1 *= -1;
		c1 *= -1;
	}
	if (s2 < 0)
	{
		a2 *= -1;
		b2 *= -1;
		c2 *= -1;
	}
	if (s3 < 0)
	{
		a3 *= -1;
		b3 *= -1;
		c3 *= -1;
	}

	/* koeficienty pro barycentricke souradnice */
	alpha = 1.0 / ABS(s2);
	beta = 1.0 / ABS(s3);
	/*gamma = 1.0 / ABS(s1);*/

	/* vyplnovani... */
	for (y = miny; y <= maxy; ++y)
	{
		/* inicilizace hranove fce v bode (minx, y) */
		e1 = a1 * minx + b1 * y + c1;
		e2 = a2 * minx + b2 * y + c2;
		e3 = a3 * minx + b3 * y + c3;

		for (x = minx; x <= maxx; ++x)
		{
			if (e1 >= 0 && e2 >= 0 && e3 >= 0)
			{
				/* interpolace pomoci barycentrickych souradnic
				e1, e2, e3 je aktualni vzdalenost bodu (x, y) od primek */
				w1 = alpha * e2;
				w2 = beta * e3;
				w3 = 1.0 - w1 - w2;



				/* interpolace zuv-souradnice */
				z = w1 * v1->z + w2 * v2->z + w3 * v3->z;
				u = w1 * t1->x + w2 * t2->x + w3 * t3->x;
				v = w1 * t1->y + w2 * t2->y + w3 * t3->y;

				//ziskani hodnoty textury
				textura = studrenTextureValue(renderer, u, v);

				//aby byly v danem rozsahu
				R = textura.red / 255.0;
				G = textura.green / 255.0;
				B = textura.blue / 255.0;

				/* interpolace barvy */
				color.red = ROUND2BYTE(w1 * col1.red + w2 * col2.red + w3 * col3.red);
				color.green = ROUND2BYTE(w1 * col1.green + w2 * col2.green + w3 * col3.green);
				color.blue = ROUND2BYTE(w1 * col1.blue + w2 * col2.blue + w3 * col3.blue);
				color.alpha = 255;

				//*hodnoty textury
				color.red = ROUND2BYTE(color.red * R);
				color.green = ROUND2BYTE(color.green * G);
				color.blue = ROUND2BYTE(color.blue * B);

				/* vykresleni bodu */
				if (z < DEPTH(pRenderer, x, y))
				{
					PIXEL(pRenderer, x, y) = color;
					DEPTH(pRenderer, x, y) = z;
				}
			}

			/* hranova fce o pixel vedle */
			e1 += a1;
			e2 += a2;
			e3 += a3;
		}
	}	
}

/******************************************************************************
 * Vykresli i-ty trojuhelnik modelu pomoci nove fce studrenDrawTriangle()
 * Pred vykreslenim aplikuje na vrcholy a normaly trojuhelniku
 * aktualne nastavene transformacni matice!
 * i - index trojuhelniku
 */

void studrenProjectTriangle(S_Renderer *pRenderer, S_Model *pModel, int i)
{
S_Coords    aa, bb, cc;             /* souradnice vrcholu po transformaci */
	S_Coords    naa, nbb, ncc;          /* normaly ve vrcholech po transformaci */
	S_Coords    nn;                     /* normala trojuhelniku po transformaci */
	int         u1, v1, u2, v2, u3, v3; /* souradnice vrcholu po projekci do roviny obrazovky */
	S_Triangle  * triangle;
	S_Coords ta, tb, tc;				/*texturovaci souradnice*/


	IZG_ASSERT(pRenderer && pModel && i >= 0 && i < trivecSize(pModel->triangles));

	/* z modelu si vytahneme trojuhelnik */
	triangle = trivecGetPtr(pModel->triangles, i);

	/* transformace vrcholu matici model */
	trTransformVertex(&aa, cvecGetPtr(pModel->vertices, triangle->v[0]));
	trTransformVertex(&bb, cvecGetPtr(pModel->vertices, triangle->v[1]));
	trTransformVertex(&cc, cvecGetPtr(pModel->vertices, triangle->v[2]));

	//texturovaci souradnice
	S_Coords *pomoc = cvecGetPtr(pModel->texcoords, triangle->v[0]);
	ta.x = pomoc->x;
	ta.y = pomoc->y;
	ta.z = pomoc->z;
	
	pomoc = cvecGetPtr(pModel->texcoords, triangle->v[1]);
	tb.x = pomoc->x;
	tb.y = pomoc->y;
	tb.z = pomoc->z;


	pomoc = cvecGetPtr(pModel->texcoords, triangle->v[2]);
	tc.x = pomoc->x;
	tc.y = pomoc->y;
	tc.z = pomoc->z;

	/* promitneme vrcholy trojuhelniku na obrazovku */
	trProjectVertex(&u1, &v1, &aa);
	trProjectVertex(&u2, &v2, &bb);
	trProjectVertex(&u3, &v3, &cc);

	/* pro osvetlovaci model transformujeme take normaly ve vrcholech */
	trTransformVector(&naa, cvecGetPtr(pModel->normals, triangle->v[0]));
	trTransformVector(&nbb, cvecGetPtr(pModel->normals, triangle->v[1]));
	trTransformVector(&ncc, cvecGetPtr(pModel->normals, triangle->v[2]));

	/* normalizace normal */
	coordsNormalize(&naa);
	coordsNormalize(&nbb);
	coordsNormalize(&ncc);

	/* transformace normaly trojuhelniku matici model */
	trTransformVector(&nn, cvecGetPtr(pModel->trinormals, triangle->n));

	/* normalizace normaly */
	coordsNormalize(&nn);

	/* je troj. privraceny ke kamere, tudiz viditelny? */
	if (!renCalcVisibility(pRenderer, &aa, &nn))
	{
		/* odvracene troj. vubec nekreslime */
		return;
	}

	/* rasterizace trojuhelniku */
        studrenDrawTriangle(pRenderer,
		&aa, &bb, &cc,
		&naa, &nbb, &ncc,
		&ta, &tb, &tc,
		u1, v1, u2, v2, u3, v3
		);


}

/******************************************************************************
 * Vrací hodnotu v aktuálnì nastavené textuøe na za
 * texturovacích souøadnicích u, v
 * Pro urèení hodnoty pou?ívá bilineární interpolaci
 * u, v - texturovací souøadnice v intervalu 0..1, který odpovídá ?íøce/vý?ce textury
 */

S_RGBA studrenTextureValue(S_StudentRenderer * pRenderer, double u, double v)
{
	if(isnan(u) || isnan(v))
		return makeColor(0, 0, 0);
	u = u * (TEXTURE_SIZE-1);
	v = v * (TEXTURE_SIZE-1);

	int x1, x2, y1, y2;
	x1 = floor(u);
	y1 = floor(v);
	x2 = ceil(u);
	y2 = ceil(v);

	double r1R, r1G, r1B, r2R, r2G, r2B;
	double pomoc1x, pomoc2x, pomoc1y, pomoc2y;

	int q11R = pRenderer->policka[x1 * TEXTURE_SIZE + y1].red;
	int q11G = pRenderer->policka[x1 * TEXTURE_SIZE + y1].green;
	int q11B = pRenderer->policka[x1 * TEXTURE_SIZE + y1].blue;

	int q12R = pRenderer->policka[x1 * TEXTURE_SIZE + y2].red;
	int q12G = pRenderer->policka[x1 * TEXTURE_SIZE + y2].green;
	int q12B = pRenderer->policka[x1 * TEXTURE_SIZE + y2].blue;

	int q21R = pRenderer->policka[x2 * TEXTURE_SIZE + y1].red;
	int q21G = pRenderer->policka[x2 * TEXTURE_SIZE + y1].green;
	int q21B = pRenderer->policka[x2 * TEXTURE_SIZE + y1].blue;

	int q22R = pRenderer->policka[x2 * TEXTURE_SIZE + y2].red;
	int q22G = pRenderer->policka[x2 * TEXTURE_SIZE + y2].green;
	int q22B = pRenderer->policka[x2 * TEXTURE_SIZE + y2].blue;

	pomoc1x = ((double)x2 - u) / ((double)x2 - x1);
	pomoc2x = (u - (double)x1) / ((double)x2 - x1);

	r1R = pomoc1x * q11R + pomoc2x * q21R;
	r1G = pomoc1x * q11G + pomoc2x * q21G;
	r1B = pomoc1x * q11B + pomoc2x * q21B;

	r2R = pomoc1x * q12R + pomoc2x * q22R;
	r2G = pomoc1x * q12G + pomoc2x * q22G;
	r2B = pomoc1x * q12B + pomoc2x * q22B;

	double pR, pG, pB;

	pomoc1y = ((double)y2 - v) / ((double)y2 - y1);
	pomoc2y = (v - (double)y1) / ((double)y2 - y1);

	pR = pomoc1y * r1R + pomoc2y * r2R;
	pG = pomoc1y * r1G + pomoc2y * r2G;
	pB = pomoc1y * r1B + pomoc2y * r2B;

	unsigned char R, G, B;
	R = ROUND2BYTE(pR);
	G = ROUND2BYTE(pG);
	B = ROUND2BYTE(pB);



    return makeColor(R, G, B);
}


/******************************************************************************
 ******************************************************************************
 * Callback funkce volana pri startu aplikace
 * Doplnte automaticke vygenerovani a prirazeni texturovacich souradnic
 * vrcholum modelu s vyuzitim mapovani na kouli
 * Muzete predpokladat, ze model je umisten v pocatku souradneho systemu
 * a posunuti neni treba resit
 */

void onInit(S_Renderer *pRenderer, S_Model *pModel)
{
	//inicializace vectoru
	vecInit(pModel->texcoords, sizeof(S_Coords));

	//delka pro cyklus
	double delka = vecSize(pModel->vertices);

	double u = 0;
	double v = 0;
	
	S_Coords vector =
	{
		.x = 0.0,
		.y = 0.0,
		.z = 0.0,

	};

	for (int i = 0; i < delka; i++)
	{
		//vytahnou souradnice
		S_Coords *souradnice = cvecGetPtr(pModel->vertices, i);

		//pozor, kralici koule
		vector.x = souradnice->x;
		vector.y = souradnice->y;
		vector.z = souradnice->z;

		//normalizace
		coordsNormalize(&vector);


		u = 0.5 + (atan2(vector.z, vector.x) / (2 * PI));
		v = 0.5 - (asin(vector.y) / PI);
		
		vector.x = u;
		vector.y = v;
		vector.z = 0;

		//a ´do vektoru
		vecPushBack(pModel->texcoords, &vector);
		

	}
}


/*****************************************************************************
 *****************************************************************************/