#include <unistd.h>
#include <fcntl.h>
short exists(char *fname){
  int fd=open(fname, O_RDONLY);
  if (fd<0)			/* error */
  return (errno==ENOENT)?-1:-2;
  /* Si no hemos salido ya, cerramos */
  close(fd);
  return 0;
}
