######################################################################
# Automatically generated by qmake (3.1) Tue Feb 6 12:39:22 2024
######################################################################

# You can make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# Please consult the documentation of the deprecated API in order to know
# how to port your code away from it.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_UP_TO=0x060000 # disables all APIs deprecated in Qt 6.0.0 and earlier

TEMPLATE = app
TARGET = homework3
INCLUDEPATH += .
QT += gui
QT += widgets
QT += core network

LIBS += -lz

# Input
SOURCES += main.cpp

# If you are using signals and slots, make sure to include the generated MOC files
# MOC_DIR is used to specify the directory where MOC will generate the files
MOC_DIR = $$OUT_PWD/moc
OBJECTS_DIR = $$OUT_PWD/obj
# Add the MOC directory to the include path
INCLUDEPATH += $$MOC_DIR
