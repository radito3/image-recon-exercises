#include <QString>
#include <QtGui>
#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QPainter>

void to_gray(QImage& image) {
    //tranformation coefficients from color to grayscale
    const double red_coef = 0.2989,
                 green_coef = 0.5870,
                 blue_coef = 0.1140;
    
    for (int row_i = 0; row_i < image.height(); row_i++) {
        QRgb* row_p = (QRgb*) image.scanLine(row_i);

        for (int col_i = 0; col_i < image.width(); col_i++) {
            auto pixel = row_p[col_i];
            int red = qRed(pixel);
            int green = qGreen(pixel);
            int blue = qBlue(pixel);

            int gray = red_coef * red + green_coef * green + blue_coef * blue;
            row_p[col_i] = qRgb(gray, gray, gray);
        }
    }
}

struct DRange {
    int min, max;
};

DRange calc_dynamic_range(const QImage& image) {
    DRange result{255, 0};

    for (int row_i = 0; row_i < image.height(); row_i++) {
        quint8* row_p = (quint8*) (image.bits() + row_i * image.bytesPerLine());

        for (int col_i = 0; col_i < image.width(); col_i++) {
            if (row_p[col_i] < result.min) {
                result.min = row_p[col_i];
            }
            if (row_p[col_i] > result.max) {
                result.max = row_p[col_i];
            }
        }
    }

    return result;
}

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    QLabel label;
    QImage image;

    const QString file_name = "lena.jpeg";

    if (image.load(file_name)) {
        QTextStream(stdout) << "Loaded: " << file_name << Qt::endl;
        
        if (image.format() == QImage::Format_RGB32) {
            to_gray(image);
            if (image.allGray()) {
                // QTextStream(stdout) << "great success!" << Qt::endl;
                QImage gray_img = image.convertToFormat(QImage::Format_Grayscale8);
                gray_img.save("gray_" + file_name);

                DRange dr = calc_dynamic_range(gray_img);
                QTextStream(stdout) << "Dynamic range: min: " << dr.min << " max: " << dr.max << Qt::endl;
                QTextStream(stdout) << "Contrast: " << dr.max - dr.min << Qt::endl;
            }
        }

        label.setPixmap(QPixmap::fromImage(image));
        label.show();
    } else {
        QTextStream(stderr) << "Failed to load: " << file_name << Qt::endl;
    }

    return app.exec();
}
