#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <omp.h>

typedef unsigned short WORD;
typedef unsigned int DWORD;
typedef long LONG; 
#define BI_RGB        0L

#pragma pack(push, 1)

typedef struct tagBITMAPFILEHEADER
{
	WORD bfType;  //specifies the file type
	DWORD bfSize;  //specifies the size in bytes of the bitmap file
	WORD bfReserved1;  //reserved; must be 0
	WORD bfReserved2;  //reserved; must be 0
	DWORD bOffBits;  //species the offset in bytes from the bitmapfileheader to the bitmap bits
}BITMAPFILEHEADER;

#pragma pack(pop)

#pragma pack(push, 1)

typedef struct tagBITMAPINFOHEADER
{
	DWORD biSize;  //specifies the number of bytes required by the struct
	LONG biWidth;  //specifies width in pixels
	LONG biHeight;  //species height in pixels
	WORD biPlanes; //specifies the number of color planes, must be 1
	WORD biBitCount; //specifies the number of bit per pixel
	DWORD biCompression;//spcifies the type of compression
	DWORD biSizeImage;  //size of image in bytes
	LONG biXPelsPerMeter;  //number of pixels per meter in x axis
	LONG biYPelsPerMeter;  //number of pixels per meter in y axis
	DWORD biClrUsed;  //number of colors used by th ebitmap
	DWORD biClrImportant;  //number of colors that are important
}BITMAPINFOHEADER;

#pragma pack(pop)

bool SaveImage(char* szPathName, void* lpBits, int w, int h) 
{
	// Create a new file for writing
	FILE* pFile = fopen(szPathName, "wb"); // wb -> w: writable b: binary, open as writable and binary
	if (pFile == NULL) 
	{
		return false;
	}

	BITMAPINFOHEADER BMIH;                         // BMP header
	BMIH.biSize = sizeof(BITMAPINFOHEADER);
	BMIH.biSizeImage = w * h * 3;
	// Create the bitmap for this OpenGL context
	BMIH.biSize = sizeof(BITMAPINFOHEADER);
	BMIH.biWidth = w;
	BMIH.biHeight = h;
	BMIH.biPlanes = 1;
	BMIH.biBitCount = 24;
	BMIH.biCompression = BI_RGB;
	BMIH.biSizeImage = w * h * 3;

	BITMAPFILEHEADER bmfh;                         // Other BMP header
	int nBitsOffset = sizeof(BITMAPFILEHEADER) + BMIH.biSize;
	LONG lImageSize = BMIH.biSizeImage;
	LONG lFileSize = nBitsOffset + lImageSize;
	bmfh.bfType = 'B' + ('M' << 8);
	bmfh.bOffBits = nBitsOffset;
	bmfh.bfSize = lFileSize;
	bmfh.bfReserved1 = bmfh.bfReserved2 = 0;

	// Write the bitmap file header               // Saving the first header to file
	DWORD nWrittenFileHeaderSize = fwrite(&bmfh, 1, sizeof(BITMAPFILEHEADER), pFile);

	// And then the bitmap info header            // Saving the second header to file
	DWORD nWrittenInfoHeaderSize = fwrite(&BMIH, 1, sizeof(BITMAPINFOHEADER), pFile);

	// Finally, write the image data itself
	//-- the data represents our drawing          // Saving the file content in lpBits to file
	DWORD nWrittenDIBDataSize = fwrite(lpBits, 1, lImageSize, pFile);
	fclose(pFile); // closing the file.

	return true;
}

#pragma region PERLIN NOISE

unsigned char permutation[] = { 151, 160, 137, 91, 90, 15,
131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23,
190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196,
135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123,
5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228,
251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107,
49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180 };

//p[i]=p[256+i]
unsigned char p[512];
void init_p()
{
	for (int i = 0; i < 256; i++)
	{
		p[i] = permutation[i];
		p[i + 256] = permutation[i];
	}
}

inline void init_p(unsigned int seed)
{
	srand(seed);
	//init half of array
	for (int i = 0; i < 256; i++)
	{
		p[i] = i;
	}
	//shuffle
	int tmp = 0, index1 = 0, index2 = 0;
	for (int i = 0; i < 1000; i++)
	{
		index1 = rand() % 256;
		index2 = rand() % 256;
		tmp = p[index1];
		p[index1] = p[index2];
		p[index2] = tmp;
	}
	//duplicate
	for (int i = 0; i < 256; i++)
	{
		p[i + 256] = p[i];
	}
}

double fade(double t)
{
	return t * t * t * (t * (t * 6 - 15) + 10);
}

double lerp(double t, double a, double b)
{
	return a + t * (b - a);
}

double grad(int hash, double x, double y, double z)
{
	int h = hash & 15;                      // CONVERT LO 4 BITS OF HASH CODE
	double u = h<8 ? x : y,                 // INTO 12 GRADIENT DIRECTIONS.
		v = h<4 ? y : h == 12 || h == 14 ? x : z;
	return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
}

double noise(double x, double y, double z) {
	int X = (int)floor(x) & 255,                  // FIND UNIT CUBE THAT
		Y = (int)floor(y) & 255,                  // CONTAINS POINT.
		Z = (int)floor(z) & 255;
	x -= floor(x);                                // FIND RELATIVE X,Y,Z
	y -= floor(y);                                // OF POINT IN CUBE.
	z -= floor(z);
	double u = fade(x),                                // COMPUTE FADE CURVES
		v = fade(y),                                // FOR EACH OF X,Y,Z.
		w = fade(z);
	int A = p[X] + Y, AA = p[A] + Z, AB = p[A + 1] + Z,      // HASH COORDINATES OF
		B = p[X + 1] + Y, BA = p[B] + Z, BB = p[B + 1] + Z;      // THE 8 CUBE CORNERS,

	return lerp(w, lerp(v, lerp(u, grad(p[AA], x, y, z),  // AND ADD
		grad(p[BA], x - 1, y, z)), // BLENDED
		lerp(u, grad(p[AB], x, y - 1, z),  // RESULTS
			grad(p[BB], x - 1, y - 1, z))),// FROM  8
		lerp(v, lerp(u, grad(p[AA + 1], x, y, z - 1),  // CORNERS
			grad(p[BA + 1], x - 1, y, z - 1)), // OF CUBE
			lerp(u, grad(p[AB + 1], x, y - 1, z - 1),
				grad(p[BB + 1], x - 1, y - 1, z - 1))));
}

#pragma endregion

#pragma region PERLIN CALLS

void SequentialPerlin(const char *fileName, BITMAPFILEHEADER &bitmapFileHeader, BITMAPINFOHEADER &bitmapInfoHeader)
{
	clock_t start = clock();
	unsigned int width, height;
	unsigned char bitmapData[3];
	FILE *filePtr; //our file pointer
				   //open filename in write binary mode
	filePtr = fopen(fileName, "wb");
	if (filePtr == NULL)
		return;
	fwrite(&bitmapFileHeader, 1, sizeof(BITMAPFILEHEADER), filePtr);
	fwrite(&bitmapInfoHeader, 1, sizeof(BITMAPINFOHEADER), filePtr);
	fseek(filePtr, bitmapFileHeader.bOffBits, SEEK_SET);
	// Visit every pixel of the image and assign a color generated with Perlin noise
	for (unsigned int i = 0; i < bitmapInfoHeader.biSizeImage; i += 3)
	{
		width = (i / 3) % bitmapInfoHeader.biWidth;
		height = i / (bitmapInfoHeader.biWidth * 3);
		double x = (double)width / ((double)bitmapInfoHeader.biWidth);
		double y = (double)height / ((double)bitmapInfoHeader.biHeight);

		// Wood like structure
		double n = 20 * noise(x, y, 0.8);
		n = n - floor(n);

		// Map the values to the [0, 255] interval, for simplicity we use tones of grey
		bitmapData[0] = floor(255 * n);
		bitmapData[1] = bitmapData[0];
		bitmapData[2] = bitmapData[0];
		fwrite(bitmapData, 1, 3, filePtr);
	}
	fclose(filePtr);
	clock_t end = clock();
	printf("czas obliczen = %f s\n", (float)(end - start) / CLOCKS_PER_SEC);
}

void ParallelPerlin(const char *fileName, BITMAPFILEHEADER &bitmapFileHeader, BITMAPINFOHEADER &bitmapInfoHeader)
{
	clock_t start = clock();
	unsigned int width, height;
	FILE *filePtr; //our file pointer
				   //open filename in write binary mode
	filePtr = fopen(fileName, "wb");
	if (filePtr == NULL)
		return;
	fwrite(&bitmapFileHeader, 1, sizeof(BITMAPFILEHEADER), filePtr);
	fwrite(&bitmapInfoHeader, 1, sizeof(BITMAPINFOHEADER), filePtr);
	fseek(filePtr, bitmapFileHeader.bOffBits, SEEK_SET);


	unsigned char* bitmapData;
	unsigned int adjustedBitmapSize;
	unsigned int shorterLoop;
	unsigned int modulo;
	unsigned int bytesToWrite;
	// (u¿ywaj¹c bariery, uwa¿aæ na ostatni¹ liniê, gdzie nie ka¿dy w¹tek mo¿e coœ obliczaæ)
	// wpisaæ masterem obliczon¹ porcjê danych do pliku (je¿eli danych bêdzie wiêcej, to master powinien uci¹æ ich czêœæ)
	#pragma omp parallel
	{
		unsigned int numberOfThreads = omp_get_num_threads();
		unsigned int threadIndex = omp_get_thread_num();
		unsigned int startIndex = threadIndex * 3;
		unsigned int loopStep = 3 * numberOfThreads;
		#pragma omp master
		{
			bitmapData = (unsigned char*)malloc(sizeof(unsigned char) * 3 * numberOfThreads);
			adjustedBitmapSize = bitmapInfoHeader.biSizeImage / 3;
			modulo = adjustedBitmapSize % numberOfThreads;
			if (modulo != 0)
			{
				adjustedBitmapSize += numberOfThreads - modulo;
			}
			adjustedBitmapSize *= 3;
			shorterLoop = loopStep - (adjustedBitmapSize - bitmapInfoHeader.biSizeImage);
		}
		#pragma omp barrier
		// Visit every pixel of the image and assign a color generated with Perlin noise
		for (unsigned int i = startIndex; i < adjustedBitmapSize; i += loopStep)
		{
			width = (i / 3) % bitmapInfoHeader.biWidth;
			height = i / (bitmapInfoHeader.biWidth * 3);
			double x = (double)width / ((double)bitmapInfoHeader.biWidth);
			double y = (double)height / ((double)bitmapInfoHeader.biHeight);

			// Wood like structure
			double n = 20 * noise(x, y, 0.8);
			n = n - floor(n);

			// Map the values to the [0, 255] interval, for simplicity we use tones of grey
			bitmapData[startIndex] = floor(255 * n);
			bitmapData[startIndex + 1] = bitmapData[startIndex];
			bitmapData[startIndex + 2] = bitmapData[startIndex];
			#pragma omp barrier
			#pragma omp master
			{
				bytesToWrite = (i + loopStep) <= bitmapInfoHeader.biSizeImage ? loopStep : shorterLoop;
				printf("bytesToWrite=%d\n", bytesToWrite);
				fwrite(bitmapData, 1, bytesToWrite, filePtr);
			}
			#pragma omp barrier
		}
	}
	free(bitmapData);
	fclose(filePtr);
	clock_t end = clock();
	printf("czas obliczen = %f s\n", (float)(end - start) / CLOCKS_PER_SEC);
}

#pragma endregion

#pragma region INITIALIZATION

void InitFileHeader(BITMAPFILEHEADER& bitmapFileHeader, unsigned int imageSize)
{
	bitmapFileHeader.bfType = 19778;
	bitmapFileHeader.bfSize = imageSize*imageSize * 3 + 122;
	bitmapFileHeader.bfReserved1 = 0;
	bitmapFileHeader.bfReserved2 = 0;
	bitmapFileHeader.bOffBits = 122;
}

void InitInfoHeader(BITMAPINFOHEADER& bitmapInfoHeader, unsigned int imageSize)
{
	bitmapInfoHeader.biSize = 108;
	bitmapInfoHeader.biWidth = imageSize;
	bitmapInfoHeader.biHeight = imageSize;
	bitmapInfoHeader.biPlanes = 1;
	bitmapInfoHeader.biBitCount = 24;
	bitmapInfoHeader.biCompression = 0;
	bitmapInfoHeader.biSizeImage = imageSize*imageSize * 3;
	bitmapInfoHeader.biXPelsPerMeter = 2835;
	bitmapInfoHeader.biYPelsPerMeter = 2835;
	bitmapInfoHeader.biClrUsed = 0;
	bitmapInfoHeader.biClrImportant = 0;
}

#pragma endregion

int main()
{
	BITMAPFILEHEADER bitmapFileHeader;
	BITMAPINFOHEADER bitmapInfoHeader;
	unsigned int seed, imageSize;
	unsigned char generationType;
	printf("Podaj ziarno do generacji wektora permutacji:\n");
	scanf("%d", &seed);
	printf("Podaj rozmiar wynikowego obrazka:\n");
	scanf("%d", &imageSize);
	init_p(seed);
	InitFileHeader(bitmapFileHeader, imageSize);
	InitInfoHeader(bitmapInfoHeader, imageSize);
	printf("Wybierz sposob obliczen: 's' - sekwencyjne, 'p' - rownolegle (OpenMP):\n");
	do
	{
		scanf("%c", &generationType);
		switch (generationType)
		{
		case 's':
			SequentialPerlin("testbitmap.bmp", bitmapFileHeader, bitmapInfoHeader);
			break;
		case 'p':
			ParallelPerlin("testbitmap.bmp", bitmapFileHeader, bitmapInfoHeader);
			break;
		}
	} while (generationType == '\n');
	system("pause");
	return 0;
}

