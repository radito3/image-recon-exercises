#include <QString>
#include <QtGui>
#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QPainter>
#include <QTextStream>

const QString file_name = "../thresh/gray_lena.jpeg";

const int INTENSITY_MIN = 0;
const int INTENSITY_MAX = 255;

void calc_histogram(const QImage& image, double histo[]) {
    for (int bean_i = 0; bean_i <= INTENSITY_MAX; bean_i++) {
        histo[bean_i] = 0;
    }
    
    int num_pixels = image.height() * image.width();
    //normalize histogram by dividing the hist. values by the num pixels
    double incr = 1.0 / num_pixels;

    for (int row_i = 0; row_i < image.height(); row_i++) {
        quint8* row_p = (quint8*) (image.bits() + row_i * image.bytesPerLine());

        for (int col_i = 0; col_i < image.width(); col_i++) {
            histo[row_p[col_i]] += incr;
        }
    }
}

void print_histogram(const double histo[]) {
    for (int bean_idx = 0; bean_idx <= INTENSITY_MAX; bean_idx++) {
        QTextStream(stdout) << '[' << bean_idx << "]: " << histo[bean_idx] << Qt::endl;
    }
}

void equalize_histogram(QImage& image, const double histo[]) {
    for (int row_i = 0; row_i < image.height(); row_i++) {
        quint8* row_p = (quint8*) (image.bits() + row_i * image.bytesPerLine());

        for (int col_i = 0; col_i < image.width(); col_i++) {
            double sum = 0;
            for (int hst_i = 0; hst_i < row_p[col_i]; hst_i++) {
                sum += histo[hst_i];
            }
            row_p[col_i] = (INTENSITY_MAX - 1) * sum;
        }
    }
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    QImage image;
    QLabel label;

    if (image.load(file_name)) {
        QTextStream(stdout) << "Image format: " << image.format() << Qt::endl;

        if (image.format() == QImage::Format_Grayscale8) {
            double histo[INTENSITY_MAX + 1];
            calc_histogram(image, histo);
            print_histogram(histo);
            equalize_histogram(image, histo);
        }
        label.setPixmap(QPixmap::fromImage(image));
        label.show();
    } else {
        QTextStream(stderr) << "Failed to load: " << file_name << Qt::endl;
    }

    return app.exec();
}
