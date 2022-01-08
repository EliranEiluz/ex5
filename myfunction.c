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

typedef struct {
    int red;
    int green;
    int blue;
    // int num;
} pixel_sum;


/*
 *  Applies kernel for pixel at (i,j)
 */

static pixel applyKernel(int dim, int i, int j, pixel *src, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

    register int ii, jj,index,index2;
    pixel_sum sum;
    pixel current_pixel,pix;
    int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
    int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
    int min_row, min_col, max_row, max_col;
    pixel loop_pixel;
    int kRow, kCol;
    register int x, y, z;
    x = y= z = 0;
    int checker = MIN(i+1, dim-1);
    int checker1 = MIN(j+1, dim-1);
    for(ii = MAX(i-1, 0); ii <= checker; ii++) {
        // compute row index in kernel
        if (ii < i) {
            kRow = 0;
        } else if (ii > i) {
            kRow = 2;
        } else {
            kRow = 1;
        }

        for(jj = MAX(j-1, 0); jj <= checker1; jj++) {


            // compute column index in kernel
            if (jj < j) {
                kCol = 0;
            } else if (jj > j) {
                kCol = 2;
            } else {
                kCol = 1;
            }

            pix = src[CALCINDEX(ii, jj, dim)];
            index = kernel[kRow][kCol];
            x = ADD(x,MULT((int) pix.red,index));
            y = ADD(y,MULT((int) pix.green,index));
            z = ADD(z,MULT((int) pix.blue,index));
        }
    }

    if (filter) {
        int toCheck;
        for(ii = MAX(i-1, 0); ii <= checker; ii++) {
            index2 = MULT(ii,dim);
            for(jj = MAX(j-1, 0); jj <= checker1; jj++) {
                // check if smaller than min or higher than max and update
                loop_pixel = src[ADD(index2,jj)];
                toCheck = ADD(ADD(((int) loop_pixel.red),((int) loop_pixel.green)),((int) loop_pixel.blue));
                if (toCheck <= min_intensity) {
                    min_intensity = toCheck;
                    min_row = ii;
                    min_col = jj;
                }
                if (toCheck > max_intensity) {
                    max_intensity = (toCheck);
                    max_row = ii;
                    max_col = jj;
                }
            }
        }
        // filter out min and max
        pix = src[CALCINDEX(min_row, min_col, dim)];
        x = ADD(x,MULT((int) pix.red,-1));
        y = ADD(y,MULT((int) pix.green,-1));
        z = ADD(z,MULT((int) pix.blue,-1));
        pix = src[CALCINDEX(max_row, max_col, dim)];
        x = ADD(x,MULT((int) pix.red,-1));
        y = ADD(y,MULT((int) pix.green,-1));
        z = ADD(z,MULT((int) pix.blue,-1));
    }
    x = DIVIDE(x, kernelScale);
    y = DIVIDE(y, kernelScale);
    z = DIVIDE(z, kernelScale);
    current_pixel.red = (unsigned char) (MIN(MAX(x, 0), 255));
    current_pixel.green = (unsigned char) (MIN(MAX(y, 0), 255));
    current_pixel.blue = (unsigned char) (MIN(MAX(z, 0), 255));
    // assign kernel's result to pixel at [i,j]
    //assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}


/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/

void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

    register int i, j, t ,w;
    int limit = SUBTRACT(dim,1);
    for (i = 1 ; i < limit; i++) {
        w = MULT(i, dim);
        for (j = 1 ; j < limit; j++) {
            //dst[calcIndex(i, j, dim)] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
            dst[ADD(w,j)] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
            //dst[ADD(w,j+1)] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
        }
    }
}


void doConvolution(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {
    int c = MULT(MULT(m,n),sizeof(pixel));
    pixel* pixelsImg = malloc(c);
    pixel* backupOrg = malloc(c);

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

    smooth(m, backupOrg, pixelsImg, kernelSize, kernel, kernelScale, filter);

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

    free(pixelsImg);
    free(backupOrg);
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

    /*
    * [1, 1, 1]
    * [1, 1, 1]
    * [1, 1, 1]
    */
    int blurKernel[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    /*
    * [-1, -1, -1]
    * [-1, 9, -1]
    * [-1, -1, -1]
    */
    int sharpKernel[3][3] = {{-1,-1,-1},{-1,9,-1},{-1,-1,-1}};

    if (flag == '1') {
        // blur image
        doConvolution(image, 3, blurKernel, 9, false);

        // write result image to file
        writeBMP(image, srcImgpName, blurRsltImgName);

        // sharpen the resulting image
        doConvolution(image, 3, sharpKernel, 1, false);

        // write result image to file
        writeBMP(image, srcImgpName, sharpRsltImgName);
    } else {
        // apply extermum filtered kernel to blur image
        doConvolution(image, 3, blurKernel, 7, true);

        // write result image to file
        writeBMP(image, srcImgpName, filteredBlurRsltImgName);

        // sharpen the resulting image
        doConvolution(image, 3, sharpKernel, 1, false);

        // write result image to file
        writeBMP(image, srcImgpName, filteredSharpRsltImgName);
    }
}