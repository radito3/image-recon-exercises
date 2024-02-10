#include <QApplication>
#include <QImage>
#include <QColor>
#include <QRandomGenerator>
#include <QMouseEvent>
#include <QPixmap>
#include <QPainter>
#include <QFileDialog>
#include <QDebug>
#include <QLabel>
#include <QVBoxLayout>
#include <QPushButton>
#include <QWidget>
#include <QRubberBand>

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <complex>

// ----------------- segmentation -----------------

// very small constant to compare floating-point numbers
const double EPS = 1.0e-14;

// gray-level intensity minimum and maximum
const int INTENS_MIN = 0;
const int INTENS_MAX = 255;

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

// accumulate normalized histogram
void calcHisto(const QImage &image, double histo[]) {
    for (int i = 0; i <= INTENS_MAX; i++) {
        histo[i] = 0;
    }

    for (int indx_row = 0; indx_row < image.height(); indx_row++) {
        quint8* ptr_row = (quint8*)(image.bits() + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++) {
            histo[ptr_row[indx_col]]++;
        }
    }

    int numb_pix = image.height() * image.width();
    for (int i = 0; i <= INTENS_MAX; i++) {
        histo[i] /= numb_pix;
    }
}

// apply threshold on a gray-scale image to calculate binary image
void thresh(QImage &image, int thr) {
    for (int indx_row = 0; indx_row < image.height(); indx_row++) {
        quint8* ptr_row = (quint8*)(image.bits() + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++) {
            ptr_row[indx_col] = (ptr_row[indx_col] < thr) ? INTENS_MIN : INTENS_MAX;
        }
    }
}

// Otsu's method
int otsu(const double histo[]) {
    // compute cumulative sums
    double p_1[INTENS_MAX + 1] = {0};
    p_1[0] = histo[0];
    for (int i = 1; i <= INTENS_MAX; i++)
    {
        p_1[i] = p_1[i - 1] + histo[i];
    }

    // cumulative mean
    double m[INTENS_MAX + 1] = {0};
    for (int i = 1; i <= INTENS_MAX; i++)
    {
        m[i] = m[i - 1] + i * histo[i];
    }

    // global mean
    double m_g = m[INTENS_MAX];

    // between-class
    double b_c[INTENS_MAX + 1] = {0};
    for (int i = 1; i <= INTENS_MAX; i++)
    {
        double div = (p_1[i] * (1 - p_1[i]));
        b_c[i] = 
            fabs(div < EPS) ? 0 :
            ((m_g * p_1[i] - m[i]) * (m_g * p_1[i] - m[i])) / div;
    }

    // find max
    double max = 0;
    int max_i = 0;
    for (int i = 0; i <= INTENS_MAX; i++)
    {
        if (b_c[i] > max)
        {
            max = b_c[i];
            max_i = i;
        }
    }

    return max_i;
}

// convert a gray-scale image to binary image
void toBinImage(QImage &image) {
    if (image.format() == QImage::Format_Grayscale8) {
        double histo[INTENS_MAX + 1];
        calcHisto(image, histo);  
        int th = otsu(histo);
        thresh(image, th);
    }
}

// ----------------- segmentation -----------------

class ImageViewer : public QWidget {
    Q_OBJECT

public:
    ImageViewer(QWidget *parent = nullptr) 
        : QWidget(parent), isImageDisplayed(false), isSelectingRect(false), isImageSegmented(false)
    {
        setupUi();
    }

protected:
    void mousePressEvent(QMouseEvent *event) override {
        if (isSelectingRect && event->button() == Qt::LeftButton) {
            selectionRect.setTopLeft(event->pos());

            origin = event->pos();
            if (!rubberBand)
                rubberBand = new QRubberBand(QRubberBand::Rectangle, this);
            rubberBand->setGeometry(QRect(origin, QSize()));
            rubberBand->show();
        }
    }

    void mouseMoveEvent(QMouseEvent *event) override {
        if (isSelectingRect) {
            selectionRect.setBottomRight(event->pos());

            rubberBand->setGeometry(QRect(origin, event->pos()).normalized());
        }
    }

    void mouseReleaseEvent(QMouseEvent *event) override {
        if (isSelectingRect && event->button() == Qt::LeftButton) {
            rubberBand->hide();

            finishCrop();
        }
    }

private slots:
    void openImage() {
        QString imagePath = QFileDialog::getOpenFileName(this, "Open Image", "", "Images (*.png *.jpg *.bmp)");

        if (!imagePath.isEmpty()) {
            displayImage(imagePath);
        }
        isImageSegmented = false;
    }

    void startCrop() {
        isSelectingRect = true;
        imageLabel->setCursor(Qt::CrossCursor);
    }

    void showSegmentedImage() {
        if (!isImageSegmented) {
            isImageSegmented = true;
            QImage binImage = transformImage(this->image);
            imageLabel->setPixmap(QPixmap::fromImage(binImage));
        }
    }

    void finishCrop() {
        isSelectingRect = false;
        isImageSegmented = false;
        imageLabel->setCursor(Qt::ArrowCursor);

        // Ensure the region of interest QRect is within the image bounds
        QRect roiRect = selectionRect.intersected(image.rect());

        QImage procssingImg = image.copy(roiRect);

        QImage binaryImg = transformImage(procssingImg);

        auto contours = detectContours(binaryImg);

        // draw contours over the original image
        QPainter painter(&image);
        painter.setPen(Qt::red);

        for (const auto& contour : contours) {
            double x = std::real(contour);
            double y = std::imag(contour);

            // Draw points on the original image within the specified QRect
            painter.drawPoint(roiRect.x() + static_cast<int>(x), roiRect.y() + static_cast<int>(y));
        }

        imageLabel->setPixmap(QPixmap::fromImage(image));
    }

private:
    std::vector<std::complex<double>> detectContours(const QImage& binaryImage) {
        std::vector<std::complex<double>> contours;

        for (int y = 0; y < binaryImage.height(); ++y) {
            for (int x = 0; x < binaryImage.width(); ++x) {
                // Check if the pixel is part of the object (black pixel in binary image)
                if (qRed(binaryImage.pixel(x, y)) == 0) {
                    // Check neighbors to determine if it is on the boundary
                    bool isBoundary = false;

                    for (int dy = -1; dy <= 1; ++dy) {
                        for (int dx = -1; dx <= 1; ++dx) {
                            int nx = x + dx;
                            int ny = y + dy;

                            // Skip pixels outside the image boundary
                            if (nx >= 0 && nx < binaryImage.width() && ny >= 0 && ny < binaryImage.height()) {
                                // Check if the neighbor is outside the object (white pixel in binary image)
                                if (qRed(binaryImage.pixel(nx, ny)) == 255) {
                                    isBoundary = true;
                                    break;
                                }
                            }
                        }

                        if (isBoundary) break;
                    }

                    // If the pixel is on the boundary, add it to the contours
                    if (isBoundary) {
                        contours.push_back(std::complex<double>(x, y));
                    }
                }
            }
        }
        return contours;
    }

    void setupUi() {
        imageLabel = new QLabel(this);
        imageLabel->setScaledContents(true);

        layout = new QVBoxLayout(this);
        layout->addWidget(imageLabel);

        openButton = new QPushButton("Open Image", this);
        layout->addWidget(openButton);

        connect(openButton, &QPushButton::clicked, this, &ImageViewer::openImage);
    }

    void displayImage(const QString &imagePath) {
        QPixmap image(imagePath);
        if (image.isNull()) {
            qDebug() << "Failed to load image:" << imagePath;
            return;
        }

        this->image = image.toImage();

        imageLabel->setPixmap(image);
        setWindowTitle("Image Viewer - " + imagePath);

        if (!isImageDisplayed) {
            isImageDisplayed = true;
            cropButton = new QPushButton("Select region", this);
            layout->addWidget(cropButton);
            segmentButton = new QPushButton("Segment image", this);
            layout->addWidget(segmentButton);

            connect(cropButton, &QPushButton::clicked, this, &ImageViewer::startCrop);
            connect(segmentButton, &QPushButton::clicked, this, &ImageViewer::showSegmentedImage);
        }
    }

    QImage transformImage(const QImage& image) {
        QImage copy = image;
        if (!image.allGray()) {
            to_gray(copy);
        }
        QImage grayImage = copy.convertToFormat(QImage::Format_Grayscale8);
        toBinImage(grayImage);
        return grayImage;
    }

    QPoint origin;
    QRubberBand* rubberBand;

    QImage image;
    QLabel *imageLabel;
    QVBoxLayout *layout;
    QPushButton *openButton;

    bool isImageDisplayed;
    QPushButton *cropButton;
    QPushButton *segmentButton;
    bool isSelectingRect;
    QRect selectionRect;
    bool isImageSegmented;
};

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    ImageViewer viewer;
    viewer.show();

    return app.exec();
}

#include "main.moc"
