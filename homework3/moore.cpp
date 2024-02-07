#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QTextStream>
#include <QString>

#include <utility> 
#include <vector>

using namespace std;

const int INTENS_MIN = 0;
const int INTENS_MAX = 255;
const int INTENS_BND = INTENS_MAX >> 1;

const QString file_name = "../bin_ligature.png";

pair<int, int> nextClockwise(const pair<int, int> &cntr, const pair<int, int> &nghb) {
    pair<int, int> result = nghb;

    int dx = nghb.first - cntr.first;
    int dy = nghb.second - cntr.second;

    if ((dx == -1 || dx == 0) && dy == -1) {
        result.first++;
    } else if (dx == 1 && (dy == -1 || dy == 0)) {
        result.second++;
    } else if ((dx == 1 || dx == 0) && dy == 1) {
        result.first--;
    } else {
        result.second--;
    }

    return result;
}

bool isBlack(const QImage &image, const pair<int, int> &pix) {
    quint8* ptr_row = (quint8*)(image.bits() 
                + pix.second * image.bytesPerLine());

    return (ptr_row[pix.first] == INTENS_MIN);
}

void trackBound(const QImage &image, vector< pair<int, int> > &bound) {
    // find uppermost leftmost point
    pair<int, int> b_0;              // uppermost leftmost point
    pair<int, int> c_0;              // west neighbor of b_0
    bool found = false;
    for (int indx_row = 0; (!found) && indx_row < image.height(); indx_row++)
    {
        quint8* ptr_row = (quint8*)(image.bits() 
                + indx_row * image.bytesPerLine());
        for (int indx_col = 0; (!found) && indx_col < image.width(); indx_col++)
        {
            if (ptr_row[indx_col] == INTENS_MIN)
            {
                found = true;
                b_0.first = indx_col;
                b_0.second = indx_row;
                c_0.first = indx_col - 1;
                c_0.second = indx_row;
            }
        }
    }

    pair<int, int> b = b_0;

    do {
        bound.push_back(b);
        pair<int, int> next = nextClockwise(b, c_0);
        while (!isBlack(image, next)) {
            c_0 = next;
            next = nextClockwise(b, c_0);
        }
        b = next;

    } while (b != b_0);
}

void printBound(const vector< pair<int, int> > &bound) {
    for (unsigned long indx_b = 0; indx_b < bound.size(); indx_b++) {
        QTextStream(stdout) << "(" << bound[indx_b].first;
        QTextStream(stdout) << ", " << bound[indx_b].second;
        QTextStream(stdout) << ") ";
    }
    QTextStream(stdout) << "\n";
}

void drawBound(QImage &image, const vector< pair<int, int> > &bound) {
    for (unsigned long indx_b = 0; indx_b < bound.size(); indx_b++) {
        quint8* ptr_row = (quint8*)(image.bits() 
                + bound[indx_b].second * image.bytesPerLine()); 
        ptr_row[bound[indx_b].first] = INTENS_BND;  
    }
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    QImage image;
    QLabel label;
    if (image.load(file_name)) {
        QTextStream(stdout) << "Image loaded: " << file_name << endl;
        QTextStream(stdout) << "Format: " << image.format() << endl;

        if (image.format() == QImage::Format_Grayscale8) {
            vector < pair<int, int> > bound;
            trackBound(image, bound);
            printBound(bound);
            drawBound(image, bound);
        } 

        label.setPixmap(QPixmap::fromImage(image));
        image.save("out.png");
    } else {
        QTextStream(stdout) << "Cannot load image: " << file_name << endl;
    }

    label.show();

    return app.exec();
}
