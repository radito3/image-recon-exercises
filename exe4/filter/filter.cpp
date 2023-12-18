#include <QString>
#include <QtGui>
#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QPainter>
#include <QTextStream>

const QString file_name = "gray_lena.jpeg";

const int INTENSITY_MIN = 0;
const int INTENSITY_MAX = 255;
const int KERNEL_SIZE = 3;

const double KERNEL_ELEM = 1.0 / (KERNEL_SIZE * KERNEL_SIZE);

double kernel[KERNEL_SIZE][KERNEL_SIZE] = {
    {KERNEL_ELEM, KERNEL_ELEM, KERNEL_ELEM},
    {KERNEL_ELEM, KERNEL_ELEM, KERNEL_ELEM},
    {KERNEL_ELEM, KERNEL_ELEM, KERNEL_ELEM}
};

//convolution (sum of multiplications) in object space:
//given a pixel map (image) of N x M,
//given a "filter" kernel with size Z x Z,
//given the root elem of the kernel be the center one (could be either the top-left, or smth. else...),
//1. multiply each pixel of the image that corresponds to the position in the kernel mask 
//2. take the sum of those multiplications
//3. write that sum in the pixel matching the position of the chosen root in the kernel mask
//4. move the filter mask along the pixel map and repeat

//example
//
//pixel map (raw image): 
//   0  1  2  3  4
//  ----------------
//0 |  |  |  |  |  |
//  ----------------
//1 |  |  |  |  |  |
//  ----------------
//2 |  |  |  |  |  |
//  ----------------
//3 |  |  |  |  |  |
//  ----------------
//
//filter mask:
// -------
// | | | |
// -------
// | |x| |
// -------
// | | | |
// -------
//
// we choose 'x' to be our kernel root
// then, we multiply the pixel intensity value (for grayscale images) of pixel (0, 0) by the kernel element
// at (0, 0) 
// iterate over the pixels that are covered by the mask and sum their value
// save that sum in the pixel that has the same coordinates with 'x' (in this case (1, 1))
// move the mask along the image and repeat the process

void filter(const QImage& in_image, QImage& out_image, double kernel[][KERNEL_SIZE]) {
    const int DK = KERNEL_SIZE / 2; //half of kernel size (non-floating point)

    for (int row_idx = 0; row_idx < in_image.height(); row_idx++) {
        quint8* out_row_ptr = (quint8*) (out_image.bits() + row_idx * out_image.bytesPerLine());
        
        for (int col_idx = 0; col_idx < out_image.width(); col_idx++) {
            double conv_value = 0;
            
            for (int ker_row_idx = 0; ker_row_idx < KERNEL_SIZE; ker_row_idx++) {
                int x = row_idx - DK + ker_row_idx;
            
                if (x >= 0 && x < in_image.height()) {
                    quint8* in_row_ptr = (quint8*) (in_image.bits() + x * in_image.bytesPerLine());
                    
                    for (int ker_col_idx = 0; ker_col_idx < KERNEL_SIZE; ker_col_idx++) {
                        int y = col_idx - DK + ker_col_idx;
                        
                        if (y >= 0 && y < in_image.width()) {
                            conv_value += in_row_ptr[y] * kernel[ker_row_idx][ker_col_idx];
                        }
                    }
                }
            }
            out_row_ptr[col_idx] = (int) conv_value;
        }
    }
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    QImage image;
    QLabel label;

    if (image.load(file_name)) {
        QTextStream(stdout) << "Image format: " << image.format() << Qt::endl;

        QImage out_image(image.width(), image.height(), QImage::Format_Grayscale8);

        if (image.format() == QImage::Format_Grayscale8) {
            filter(image, out_image, kernel);
        }

        label.setPixmap(QPixmap::fromImage(out_image));
        label.show();
    } else {
        QTextStream(stderr) << "Failed to load: " << file_name << Qt::endl;
    }

    return app.exec();
}
