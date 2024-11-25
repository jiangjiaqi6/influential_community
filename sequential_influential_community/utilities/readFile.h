// #ifndef READFILE
// #define READFILE
// #include "Defines.h"
// constexpr uint32_t INF = 0x7fffffff;

// class readFile{
// private:
//     long len;
//     int fd;
//     char *addr;
//     char *bptr,*eptr;
//     bool isFinish;

// public:
//     readFile(const char* path){
//         len = 0, fd = -1;
//         isFinish = false;
//         addr = NULL;
//         struct stat statbuf;
//         fd = open(path,O_RDONLY);
//         assert(fd >= 0 );
//         int ret = fstat(fd, &statbuf);
//         assert(ret >= 0);
//         len = statbuf.st_size;
//         addr = (char*)mmap(0, len, PROT_READ, MAP_SHARED, fd, 0);
//         assert (addr != (void *)-1);
//         bptr = addr;
//         eptr = addr + len;
//     }
//     ~readFile(){
//         if(fd != -1)
//         {
//             close(fd);
//             munmap(addr, len);
//         }
//     }

//     unsigned int getUInt(){
//         uint32_t ret = 0;
//         bool readF = false;
        
//         while (bptr < eptr)
//         {
//             char ch = *bptr;
//             if (ch >= '0' && ch <= '9')
//             {
//                 ret = ret * 10 + (ch - '0');
//                 readF = true;
//             }
//             else
//             {
//                 if (readF)
//                     return ret;
//             }
//             bptr++;
//         }
//         return ret;
//     }

//     float getFloat(){
//         float right = 0, left = 0;
//         char ch;
//         while (bptr < eptr && (*bptr < '0' || *bptr > '9'))
//             bptr++;
//         while(bptr < eptr && (*bptr >= '0' && *bptr <= '9')){
//             left = left * 10 + (*bptr - '0');
//             bptr++;
//         }
//         if(*bptr == '.'){
//             for(ch = *bptr; bptr < eptr && (ch < '0' || ch > '9'); ) 
//                 bptr++,ch = *bptr;
//             float base = 10;
//             while(bptr < eptr && (*bptr >= '0' && *bptr <= '9')){
//                 left += float(*bptr - '0')/base;
//                 bptr++;
//                 base *= 10;
//             }
//         }        
//         return left;
//     }

//     bool empty(){
//         if(*bptr == '%'){
//             while(*bptr != '\n')
//                 bptr++;
//             bptr++;
//             if(*bptr == '%'){
//                 while(*bptr != '\n')
//                     bptr++;
//             }
//         }
//         while(bptr < eptr && (*bptr < '0' || *bptr > '9')) bptr++;
//         if(bptr == eptr) isFinish = true;
//         return isFinish;
//     }

//     long getLen()
//     {
//         return len;
//     }

//     char *getAddr()
//     {
//         return addr;
//     }

// };
// #endif