#include <stdbool.h>
#define MAX(a,b) ((a) > (b) ? (a):(b))
#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MULT(a,b) ((a) * (b))
#define SUBTRACT(a,b) ((a) - (b))
#define ADD(a,b) ((a) + (b))
#define DIVIDE(a,b) ((a) / (b))
#define CALCINDEX(a,b,c) (((a) * (c)) + (b))

typedef struct {
    unsigned char red;
    unsigned char green;
    unsigned char blue;
} pixel;


/*
 *  Applies kernel for pixel at (i,j)
 */


void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

    int c = MULT(MULT(m,n),sizeof(pixel));
    pixel* pixelsImg = malloc(c);
    pixel* backupOrg = malloc(c);
    register int i, j, w1,w2,w3;
    int limit = SUBTRACT(m,1);
    register int sum1, sum2, sum3, sum4;
    register int row, col;
    for (row = 0 ; row < m ; row++) {
        sum1 = MULT(row,n);
        sum2 = MULT(sum1,3);
        for (col = 0 ; col < n ; col++) {
            sum4 = ADD(MULT(3,col),sum2);
            sum3 = ADD(sum1,col);
            backupOrg[sum3].red = pixelsImg[sum3].red = image->data[sum4];
            backupOrg[sum3].green = pixelsImg[sum3].green = image->data[ADD(sum4, 1)];
            backupOrg[sum3].blue = pixelsImg[sum3].blue = image->data[ADD(sum4, 2)];
        }
    }

    if (flag == '1') {
        // blur image
        //doConvolution(image, 3, blurKernel, 9, false);
        for (i = 1 ; i < limit; i++) {
            w1 = MULT(i-1,m);
            w2 = MULT(i, m);
            w3 = MULT(i+1,m);
            for (j = 1 ; j < limit; j++) {
                pixelsImg[ADD(w2,j)].red = DIVIDE((backupOrg[ADD(w1,j-1)].red + backupOrg[ADD(w1,j)].red + backupOrg[ADD(w1,j+1)].red
                                                   + backupOrg[ADD(w2,j-1)].red + backupOrg[ADD(w2,j)].red + backupOrg[ADD(w2,j+1)].red +
                                                   backupOrg[ADD(w3,j-1)].red + backupOrg[ADD(w3,j)].red + backupOrg[ADD(w3,j+1)].red),9);

                pixelsImg[ADD(w2,j)].green = DIVIDE((backupOrg[ADD(w1,j-1)].green + backupOrg[ADD(w1,j)].green + backupOrg[ADD(w1,j+1)].green +
                                                   backupOrg[ADD(w2,j-1)].green + backupOrg[ADD(w2,j)].green + backupOrg[ADD(w2,j+1)].green +
                                                   backupOrg[ADD(w3,j-1)].green + backupOrg[ADD(w3,j)].green + backupOrg[ADD(w3,j+1)].green),9);

                pixelsImg[ADD(w2,j)].blue = DIVIDE((backupOrg[ADD(w1,j-1)].blue + backupOrg[ADD(w1,j)].blue + backupOrg[ADD(w1,j+1)].blue
                                                   + backupOrg[ADD(w2,j-1)].blue + backupOrg[ADD(w2,j)].blue + backupOrg[ADD(w2,j+1)].blue +
                                                   backupOrg[ADD(w3,j-1)].blue + backupOrg[ADD(w3,j)].blue + backupOrg[ADD(w3,j+1)].blue),9);
            }
        }
        for (row = 0 ; row < m ; row++) {
            sum1 = MULT(row,n);
            sum2 = MULT(sum1,3);
            for (col = 0 ; col < n ; col++) {
                sum4 = ADD(MULT(3,col),sum2);
                sum3 = ADD(sum1,col);
                backupOrg[sum3].red = image->data[sum4] = pixelsImg[sum3].red;
                backupOrg[sum3].green = image->data[ADD(sum4,1)] = pixelsImg[sum3].green;
                backupOrg[sum3].blue = image->data[ADD(sum4,2)] = pixelsImg[sum3].blue;
            }
        }
        // write result image to file
        writeBMP(image, srcImgpName, blurRsltImgName);

    } else {
        int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
        int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
        // apply extermum filtered kernel to blur image
        //doConvolution(image, 3, blurKernel, 7, true);
        int index2,toCheck,min,max;
        pixel loop_pixel;
        for (i = 1 ; i < limit; i++) {
            w1 = MULT(i-1,m);
            w2 = MULT(i, m);
            w3 = MULT(i+1,m);
            for (j = 1 ; j < limit; j++) {
                for (int ii = i-1; ii <= i+1; ii++) {
                    index2 = MULT(ii, m);
                    for (int jj = j-1; jj <= j+1; jj++) {
                        // check if smaller than min or higher than max and update
                        loop_pixel = backupOrg[ADD(index2, jj)];
                        toCheck = ADD(ADD((loop_pixel.red), (loop_pixel.green)), (loop_pixel.blue));
                        if (toCheck <= min_intensity) {
                            min_intensity = toCheck;
                            min = index2 +jj;
                        }
                        if (toCheck > max_intensity) {
                            max_intensity = (toCheck);
                            max = index2 + jj;
                        }
                    }
                }
                //dst[calcIndex(i, j, dim)] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
                pixelsImg[ADD(w2,j)].red = DIVIDE((backupOrg[ADD(w1,j-1)].red + backupOrg[ADD(w1,j)].red + backupOrg[ADD(w1,j+1)].red
                                                   + backupOrg[ADD(w2,j-1)].red + backupOrg[ADD(w2,j)].red + backupOrg[ADD(w2,j+1)].red +
                                                   backupOrg[ADD(w3,j-1)].red + backupOrg[ADD(w3,j)].red + backupOrg[ADD(w3,j+1)].red - backupOrg[min].red -backupOrg[max].red),7);
                pixelsImg[ADD(w2,j)].green = DIVIDE((backupOrg[ADD(w1,j-1)].green + backupOrg[ADD(w1,j)].green + backupOrg[ADD(w1,j+1)].green
                                                     + backupOrg[ADD(w2,j-1)].green + backupOrg[ADD(w2,j)].green + backupOrg[ADD(w2,j+1)].green +
                                                     backupOrg[ADD(w3,j-1)].green + backupOrg[ADD(w3,j)].green + backupOrg[ADD(w3,j+1)].green - backupOrg[min].green -backupOrg[max].green),7);
                pixelsImg[ADD(w2,j)].blue = DIVIDE((backupOrg[ADD(w1,j-1)].blue + backupOrg[ADD(w1,j)].blue + backupOrg[ADD(w1,j+1)].blue
                                                    + backupOrg[ADD(w2,j-1)].blue + backupOrg[ADD(w2,j)].blue + backupOrg[ADD(w2,j+1)].blue +
                                                    backupOrg[ADD(w3,j-1)].blue + backupOrg[ADD(w3,j)].blue + backupOrg[ADD(w3,j+1)].blue - backupOrg[min].blue - backupOrg[max].blue),7);

            }
        }
        for (row = 0 ; row < m ; row++) {
            sum1 = MULT(row,n);
            sum2 = MULT(sum1,3);
            for (col = 0 ; col < n ; col++) {
                sum4 = ADD(MULT(3,col),sum2);
                sum3 = ADD(sum1,col);
                backupOrg[sum3].red = image->data[sum4] = pixelsImg[sum3].red;
                backupOrg[sum3].green = image->data[ADD(sum4,1)] = pixelsImg[sum3].green;
                backupOrg[sum3].blue = image->data[ADD(sum4,2)] = pixelsImg[sum3].blue;
            }
        }
        // write result image to file
        writeBMP(image, srcImgpName, filteredBlurRsltImgName);
    }

    //end of else statement

    // sharpen the resulting image
    //doConvolution(image, 3, sharpKernel, 1, false);
    for (i = 1 ; i < limit; i++) {
        w1 = MULT(i-1,m);
        w2 = MULT(i, m);
        w3 = MULT(i+1,m);
        for (j = 1 ; j < limit; j++) {
            //dst[calcIndex(i, j, dim)] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
            pixelsImg[ADD(w2,j)].red = MIN(MAX(-1*backupOrg[ADD(w1,j-1)].red - backupOrg[ADD(w1,j)].red - backupOrg[ADD(w1,j+1)].red
                                               - backupOrg[ADD(w2,j-1)].red + 9*backupOrg[ADD(w2,j)].red - backupOrg[ADD(w2,j+1)].red -
                                               backupOrg[ADD(w3,j-1)].red - backupOrg[ADD(w3,j)].red - backupOrg[ADD(w3,j+1)].red,0),255);
            pixelsImg[ADD(w2,j)].green = MIN(MAX(-1*backupOrg[ADD(w1,j-1)].green - backupOrg[ADD(w1,j)].green - backupOrg[ADD(w1,j+1)].green
                                                 - backupOrg[ADD(w2,j-1)].green + 9*backupOrg[ADD(w2,j)].green - backupOrg[ADD(w2,j+1)].green -
                                                 backupOrg[ADD(w3,j-1)].green - backupOrg[ADD(w3,j)].green - backupOrg[ADD(w3,j+1)].green,0),255);
            pixelsImg[ADD(w2,j)].blue = MIN(MAX(-1*backupOrg[ADD(w1,j-1)].blue - backupOrg[ADD(w1,j)].blue - backupOrg[ADD(w1,j+1)].blue
                                                - backupOrg[ADD(w2,j-1)].blue + 9*backupOrg[ADD(w2,j)].blue - backupOrg[ADD(w2,j+1)].blue -
                                                backupOrg[ADD(w3,j-1)].blue - backupOrg[ADD(w3,j)].blue - backupOrg[ADD(w3,j+1)].blue,0),255);
        }
    }
    for (row = 0 ; row < m ; row++) {
        sum1 = MULT(row,n);
        sum2 = MULT(sum1,3);
        for (col = 0 ; col < n ; col++) {
            sum4 = ADD(MULT(3,col),sum2);
            sum3 = ADD(sum1,col);
            image->data[sum4] = pixelsImg[sum3].red;
            image->data[ADD(sum4,1)] = pixelsImg[sum3].green;
            image->data[ADD(sum4,2)] = pixelsImg[sum3].blue;
        }
    }

    if(flag == '1')
        writeBMP(image, srcImgpName, sharpRsltImgName);
    else
        writeBMP(image, srcImgpName, filteredSharpRsltImgName);
    free(pixelsImg);
    free(backupOrg);
}