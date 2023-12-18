#include <QString>
#include <QtGui>
#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QPainter>
#include <QTextStream>

#include <cmath>

const QString file_name = "gray_lena.jpeg";

const int INTENSITY_MIN = 0;
const int INTENSITY_MAX = 255;
const int DEFAULT_THRESHOLD = 100;

void format_by_threshold(QImage& image, int threshold) {
    for (int row_i = 0; row_i < image.height(); row_i++) {
        quint8* row_p = (quint8*) (image.bits() + row_i * image.bytesPerLine());

        for (int col_i = 0; col_i < image.width(); col_i++) {
            row_p[col_i] = row_p[col_i] < threshold ? INTENSITY_MIN : INTENSITY_MAX;
        }
    }
}

//negative image
void invert_intensity(QImage& image) {
    for (int row_i = 0; row_i < image.height(); row_i++) {
        quint8* row_p = (quint8*) (image.bits() + row_i * image.bytesPerLine());

        for (int col_i = 0; col_i < image.width(); col_i++) {
            row_p[col_i] = INTENSITY_MAX - row_p[col_i];
        }
    }
}

void log_transform(QImage& image) {
    for (int row_i = 0; row_i < image.height(); row_i++) {
        quint8* row_p = (quint8*) (image.bits() + row_i * image.bytesPerLine());

        for (int col_i = 0; col_i < image.width(); col_i++) {
            row_p[col_i] = DEFAULT_THRESHOLD * std::log(1 + row_p[col_i]);
        }
    }
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    QImage image;
    QLabel label;

    if (image.load(file_name)) {
        QTextStream(stdout) << "Image format: " << image.format() << Qt::endl;
        // QTextStream(stdout) << "Threshold: ";
        // int threshold;
        // QTextStream(stdin) >> threshold;

        if (image.format() == QImage::Format_Grayscale8) {
            // format_by_threshold(image, threshold);
            // invert_intensity(image);
            log_transform(image);
        }
        label.setPixmap(QPixmap::fromImage(image));
        label.show();
    } else {
        QTextStream(stderr) << "Failed to load: " << file_name << Qt::endl;
    }

    return app.exec();
}
