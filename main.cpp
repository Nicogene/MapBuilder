#include "MapBuilder.h"
#include <iostream>

int main()
{
    MapBuilder m(10,true,CVSBA);
    while(true){

        if(!m.process())
            break;
    }
}
