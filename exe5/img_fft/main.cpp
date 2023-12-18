#include <QString>
#include <QtGui>
#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QPainter>
#include <QTextStream>

#include <vector>
#include <complex>
#include <cmath>

using namespace std;

const QString file_name = "gray_lena.jpeg";

const int INTENSITY_MIN = 0;
const int INTENSITY_MAX = 255;

void dft2d(const QImage& image, vector<vector<complex<double>>>& trs) {
    //dimensions
    int m = image.height();
    int n = image.width();

    //default-contstructed elements are inserted
    trs.resize(m);
    for (int row = 0; row < m; ++row) {
        trs[row].resize(n);
    }

    for (int trs_row = 0; trs_row < m; trs_row++) {
        for (int trs_col = 0; trs_col < n; trs_col++) {
            complex<double> sum(0, 0);
            for (int img_row = 0; img_row < m; img_row++) {
                quint8* row_ptr = (quint8*) (image.bits() + img_row * image.bytesPerLine());
                for (int img_col = 0; img_col < n; img_col++) {
                    double expnt = -2 * M_PI * (
                        1.0 * trs_row * img_row / m 
                        + 1.0 * trs_col * img_col / n);
                    double centered = row_ptr[img_col] * pow(-1.0, img_row + img_col);
                    sum += centered * exp(complex<double>(0, expnt));
                }
            }
            trs[trs_row][trs_col] = sum;
        }
    }
}

void calc_spectrum(/*vector*/) {
    //TODO: implement
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    QImage image;
    QLabel label;

    if (image.load(file_name)) {
        QTextStream(stdout) << "Image format: " << image.format() << Qt::endl;

        QImage spectrum(image.width(), image.height(), QImage::Format_Grayscale8);

        if (image.format() == QImage::Format_Grayscale8) {
            vector<vector<complex<double>>> transformed;
            dft2d(image, transformed);
            // calc_spectrum(transformed);
        }

        label.setPixmap(QPixmap::fromImage(spectrum));
        label.show();
    } else {
        QTextStream(stderr) << "Failed to load: " << file_name << Qt::endl;
    }

    return app.exec();
}
