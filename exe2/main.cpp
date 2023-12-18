#include <QString>
#include <QtGui>
#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QPainter>

struct PCoord {
    int row, col;
};

PCoord center(QImage& image) {
    const int rect_size = 20;
    QPainter painter(&image);
    painter.setBrush(Qt::NoBrush);
    painter.setPen(Qt::green);
    PCoord center;
    center.row = image.height() / 2;
    center.col = image.width() / 2;
    painter.drawRect(center.row, center.col, rect_size, rect_size);
    painter.end();
    return center;
}

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    QLabel label;
    QImage image;

    const QString file_name = "lena.jpeg";

    if (image.load(file_name)) {
        QTextStream(stdout) << "Loaded: " << file_name << Qt::endl;
        PCoord cnr = center(image);
        QTextStream(stdout) << "row: " << cnr.row << " col: " << cnr.col << Qt::endl;
        label.setPixmap(QPixmap::fromImage(image));
        label.show();
    } else {
        QTextStream(stderr) << "Failed to load: " << file_name << Qt::endl;
    }

    return app.exec();
}
