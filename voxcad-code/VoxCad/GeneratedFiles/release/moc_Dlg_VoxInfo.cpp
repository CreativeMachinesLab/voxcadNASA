/****************************************************************************
** Meta object code from reading C++ file 'Dlg_VoxInfo.h'
**
** Created: Tue Sep 24 13:43:58 2013
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../VoxCad/Dlg_VoxInfo.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'Dlg_VoxInfo.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Dlg_VoxInfo[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: signature, parameters, type, tag, flags
      22,   13,   12,   12, 0x05,
      45,   40,   12,   12, 0x05,
      82,   71,   12,   12, 0x05,

 // slots: signature, parameters, type, tag, flags
     113,   12,   12,   12, 0x0a,
     124,   12,   12,   12, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_Dlg_VoxInfo[] = {
    "Dlg_VoxInfo\0\0CurIndex\0GetCurIndex(int*)\0"
    "Info\0GetDMInfoString(QString*)\0"
    "Index,Info\0GetVoxInfoString(int,QString*)\0"
    "UpdateUI()\0UpdateText()\0"
};

void Dlg_VoxInfo::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        Dlg_VoxInfo *_t = static_cast<Dlg_VoxInfo *>(_o);
        switch (_id) {
        case 0: _t->GetCurIndex((*reinterpret_cast< int*(*)>(_a[1]))); break;
        case 1: _t->GetDMInfoString((*reinterpret_cast< QString*(*)>(_a[1]))); break;
        case 2: _t->GetVoxInfoString((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< QString*(*)>(_a[2]))); break;
        case 3: _t->UpdateUI(); break;
        case 4: _t->UpdateText(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData Dlg_VoxInfo::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Dlg_VoxInfo::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_Dlg_VoxInfo,
      qt_meta_data_Dlg_VoxInfo, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Dlg_VoxInfo::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Dlg_VoxInfo::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Dlg_VoxInfo::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Dlg_VoxInfo))
        return static_cast<void*>(const_cast< Dlg_VoxInfo*>(this));
    return QWidget::qt_metacast(_clname);
}

int Dlg_VoxInfo::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    }
    return _id;
}

// SIGNAL 0
void Dlg_VoxInfo::GetCurIndex(int * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void Dlg_VoxInfo::GetDMInfoString(QString * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void Dlg_VoxInfo::GetVoxInfoString(int _t1, QString * _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}
QT_END_MOC_NAMESPACE
