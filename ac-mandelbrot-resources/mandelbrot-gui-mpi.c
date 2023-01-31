//mpirun -np 2 --hostfile localhost.OPENMPI ./mandelbrot-gui-mpi.exe

#include <sys/time.h> 
struct timeval start, end;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include <mpi.h>

////////////////////////////////////////////////////////////////////////
//user defined datatypes
typedef struct {unsigned char r, g, b;} rgb_t;

////////////////////////////////////////////////////////////////////////
//global variables
//(do not change from the defaults except if solicited in the work text)

rgb_t **GLOBAL_tex = 0;
int GLOBAL_gwin;
GLuint GLOBAL_texture;
int GLOBAL_width, GLOBAL_height;
int GLOBAL_tex_w, GLOBAL_tex_h;
double GLOBAL_scale = 1./256;
double GLOBAL_cx = -.6, GLOBAL_cy = 0;
int GLOBAL_color_rotate = 0;
int GLOBAL_saturation = 1;
int GLOBAL_invert = 0;

// GLOBAL_zoomin contains the path (sequence of <x,y> mouse coordinates) to zoom in the fractal)
#define GLOBAL_zoomin_num_pairs 36
int GLOBAL_zoomin[GLOBAL_zoomin_num_pairs]={538,237,491,369,522,383,492,372,504,369,510,385,491,353,472,329,516,392,518,393,478,327,506,400,501,320,487,361,501,363,501,363,486,363,486,363};

int GLOBAL_window_width=1024;
int GLOBAL_window_height=768;
int GLOBAL_refresh=1;
int GLOBAL_max_iter = 4096; // 4096 ou 256
int GLOBAL_tex_size=0;

// Tornamos globais
int numtasks, rank;

////////////////////////////////////////////////////////////////////////
//function prototypes
void render();
void keypress(unsigned char key, int x, int y);
void mouseclick(int button, int state, int x, int y);
void resize(int w, int h);
void init_gfx(int *c, char **v, int rank);

void alloc_tex();
void set_texture();

void hsv_to_rgb(int hue, int min, int max, rgb_t *p);
void calc_mandel();

void print_menu();
void screen_dump();

////////////////////////////////////////////////////////////////////////
// function implementations

void print_menu()
{
	printf("\n\nkeys:\n\t"
	"q: quit	\n\t"
	"ESC: reset to initial frame\n\t"
	"r: color rotation\n\t"
	"c: monochrome\n\t"
	"s: screen dump\n\t"
	"<, >: decrease/increase max iteration\n\t"
	"I: max iteration=4096\n\t"
	"i: max iteration=128\n\t"
	"mouse buttons to zoom\n\t"
	"z: automatic zoom in (one step)\n\t"
	"Z: automatic zoom in (all steps)\n");
}

//////////////////////////////////////////////////////////////////////// 
void screen_dump()
{
	// TERMINAR DE CONTAR O TEMPO E PRINTAR
	gettimeofday(&end,NULL);
	int elapsed = ((end.tv_sec - start.tv_sec) * 1000000) + (end.tv_usec - start.tv_usec);
	printf("\nTime elapsed: %d micro seconds\n",elapsed);

	static int dump=1;
	char fn[100];
	int i;
	sprintf(fn, "screen%03d-mpi.ppm", dump++);
	FILE *fp = fopen(fn, "w");
	fprintf(fp, "P6\n%d %d\n255\n", GLOBAL_width, GLOBAL_height);
	for (i = GLOBAL_height - 1; i >= 0; i--)
		fwrite(GLOBAL_tex[i], 1, GLOBAL_width * 3, fp);
	fclose(fp);
	printf("%s written\n", fn);
}

////////////////////////////////////////////////////////////////////////
void hsv_to_rgb(int hue, int min, int max, rgb_t *p)
{
	if (min == max) max = min + 1;
	if (GLOBAL_invert) hue = max - (hue - min);
	if (!GLOBAL_saturation) {
		p->r = p->g = p->b = 255 * (max - hue) / (max - min);
		return;
	}
	double h = fmod(GLOBAL_color_rotate + 1e-4 + 4.0 * (hue - min) / (max - min), 6);
#	define VAL 255
	double c = VAL * GLOBAL_saturation;
	double X = c * (1 - fabs(fmod(h, 2) - 1));
 
	p->r = p->g = p->b = 0;
 
	switch((int)h) {
	case 0: p->r = c; p->g = X; return;
	case 1:	p->r = X; p->g = c; return;
	case 2: p->g = c; p->b = X; return;
	case 3: p->g = X; p->b = c; return;
	case 4: p->r = X; p->b = c; return;
	default:p->r = c; p->b = X;
	}
}

int bc_array_int[12];
double bc_array_double[3];

void into_arrays(){
	bc_array_int[0] = GLOBAL_window_width;
	bc_array_int[1] = GLOBAL_window_height;
	bc_array_int[2] = GLOBAL_tex_h;
	bc_array_int[3] = GLOBAL_tex_w;
	bc_array_int[4] = GLOBAL_refresh;
	bc_array_int[5] = GLOBAL_width;
	bc_array_int[6] = GLOBAL_height;
	bc_array_int[7] = GLOBAL_max_iter;
	bc_array_int[8] = GLOBAL_invert;
	bc_array_int[9] = GLOBAL_saturation;
	bc_array_int[10] = GLOBAL_color_rotate;
	bc_array_int[11] = GLOBAL_tex_size;

	bc_array_double[0] = GLOBAL_scale;
	bc_array_double[1] = GLOBAL_cy;
	bc_array_double[2] = GLOBAL_cx;
}

void out_arrays(){
	GLOBAL_window_width = bc_array_int[0];
	GLOBAL_window_height = bc_array_int[1];
	GLOBAL_tex_h = bc_array_int[2];
	GLOBAL_tex_w = bc_array_int[3];
	GLOBAL_refresh = bc_array_int[4];
	GLOBAL_width = bc_array_int[5];
	GLOBAL_height = bc_array_int[6];
	GLOBAL_max_iter = bc_array_int[7];
	GLOBAL_invert = bc_array_int[8];
	GLOBAL_saturation = bc_array_int[9];
	GLOBAL_color_rotate = bc_array_int[10];
	GLOBAL_tex_size = bc_array_int[11];

	GLOBAL_scale = bc_array_double[0];
	GLOBAL_cy = bc_array_double[1];
	GLOBAL_cx = bc_array_double[2];
}

////////////////////////////////////////////////////////////////////////
void calc_mandel() 
{
	int i;

	if(rank == 0)
	{
		into_arrays(); // Passa todas as variaveis para dois arrays
	}

	MPI_Bcast(&bc_array_int, 12, MPI_INT, 0, MPI_COMM_WORLD); // Manda todos os ints
	MPI_Bcast(&bc_array_double, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Manda todos os doubles

	if(rank != 0) 
	{	
		out_arrays(); // Passa todos os elementos dos arrays para as variaveis
		GLOBAL_tex = realloc(GLOBAL_tex, GLOBAL_tex_size); // Alocar tamanho pro global_tex
	}
	MPI_Bcast(GLOBAL_tex, GLOBAL_tex_size, MPI_BYTE, 0, MPI_COMM_WORLD); // Manda todos os bytes do GLOBAL_tex


	if(rank != 0) // Recalcular enderecos das linhas
	{
		for (GLOBAL_tex[0] = (rgb_t *)(GLOBAL_tex + GLOBAL_tex_h), i = 1; i < GLOBAL_tex_h; i++)
			GLOBAL_tex[i] = GLOBAL_tex[i - 1] + GLOBAL_tex_w;
	}


	int j, iter, min, max;
	rgb_t *px; // nosso ponteiro local para acessar o pixel
	double x, y, zx, zy, zx2, zy2;
	min = GLOBAL_max_iter; max = 0;

	// Quantas cada task vai fazer (só pode ser par)
	int linhas_cada = GLOBAL_height/numtasks;
	int lim_inf = rank*linhas_cada;
	int lim_sup = lim_inf + linhas_cada;

	// Calcula o px, min, max
	for (i = lim_inf; i < lim_sup; i++) {
		// q é uma sequencia de variaveis rbg (3 chars) 
		px = GLOBAL_tex[i]; 

		y = (i - GLOBAL_height/2) * GLOBAL_scale + GLOBAL_cy;
		for (j = 0; j  < GLOBAL_width; j++, px++) {
			x = (j - GLOBAL_width/2) * GLOBAL_scale + GLOBAL_cx;
			iter = 0;
 
			zx = hypot(x - .25, y);
			if (x < zx - 2 * zx * zx + .25) iter = GLOBAL_max_iter;
			if ((x + 1)*(x + 1) + y * y < 1/16) iter = GLOBAL_max_iter;
 
			zx = zy = zx2 = zy2 = 0;
			for (; iter < GLOBAL_max_iter && zx2 + zy2 < 4; iter++) {
				zy = 2 * zx * zy + y;
				zx = zx2 - zy2 + x;
				zx2 = zx * zx;
				zy2 = zy * zy;
			}
			if (iter < min) min = iter;
			if (iter > max) max = iter;

			*(unsigned short *)px = iter; // Atraves de ponteiro ta alterando o valor armazenado naquele endereço
		}
	}

	// Passa min, max e cada px para a tela
	for (i = lim_inf; i < lim_sup; i++)
		for (j = 0, px = GLOBAL_tex[i]; j  < GLOBAL_width; j++, px++)
			hsv_to_rgb(*(unsigned short*)px, min, max, px); // Passa o valor atraves de ponteiro e tambem o endereco 

	// Buffer para envio
	rgb_t *sendbuf = 0;
	
	// Tamanho do buffer a ser enviado	
	int SEND_BUFFER_SIZE = linhas_cada * GLOBAL_tex_w * sizeof(rgb_t); 
	// Espaco alocado para o buffer de envio
	sendbuf = realloc(sendbuf, SEND_BUFFER_SIZE);


	int aux = 0;
	// Preenche o buffer de envio com todos os pixeis que a task calculou
	for (i = lim_inf; i < lim_sup; i++)
	{
		for (j = 0, px = GLOBAL_tex[i]; j < GLOBAL_width; j++, px++)
		{
			sendbuf[aux] = *px;
			aux++;
		}
	}

	// Envia cada pedaco do GLOBAL_tex calculado para a task 0
	MPI_Gather(sendbuf, SEND_BUFFER_SIZE, MPI_BYTE, *GLOBAL_tex, SEND_BUFFER_SIZE, MPI_BYTE, 0, MPI_COMM_WORLD);

}

////////////////////////////////////////////////////////////////////////
void alloc_tex()
{
	int i, ow = GLOBAL_tex_w, oh = GLOBAL_tex_h;
	       //backup current texture dimensions

    // if necessary, adjust texture dimensions to image dimensions 
	for (GLOBAL_tex_w = 1; GLOBAL_tex_w < GLOBAL_width;  GLOBAL_tex_w <<= 1);
	for (GLOBAL_tex_h = 1; GLOBAL_tex_h < GLOBAL_height; GLOBAL_tex_h <<= 1);
 
    // if texture dimensions were changed, realloc the texture memory block (GLOBAL_tex)
	if (GLOBAL_tex_h != oh || GLOBAL_tex_w != ow) {
	    GLOBAL_tex_size = GLOBAL_tex_h * sizeof(rgb_t*) + GLOBAL_tex_h * GLOBAL_tex_w * 3;
		GLOBAL_tex = realloc(GLOBAL_tex, GLOBAL_tex_size);
		// a texture has GLOBAL_tex_h pointers to the begining of each line,
		// followed by the GLOBAL_tex_h lines of the image (bottom to top);
        // each line is sequence of GLOBAL_tex_w*3 bytes (each pixel RGB values)
    }
  
    // update pointers in the beggining of the texture
	for (GLOBAL_tex[0] = (rgb_t *)(GLOBAL_tex + GLOBAL_tex_h), i = 1; i < GLOBAL_tex_h; i++)
		GLOBAL_tex[i] = GLOBAL_tex[i - 1] + GLOBAL_tex_w;
		// uses rgb_t* arithmetic pointer, where each unit corresponds to 3 bytes
}

////////////////////////////////////////////////////////////////////////
void set_texture()
{
	if (GLOBAL_refresh)
		alloc_tex();

	calc_mandel();
 
	if (GLOBAL_refresh) {
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, GLOBAL_texture);
		glTexImage2D(GL_TEXTURE_2D, 0, 3, GLOBAL_tex_w, GLOBAL_tex_h,
			0, GL_RGB, GL_UNSIGNED_BYTE, GLOBAL_tex[0]);
 
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		render();
	}
}

////////////////////////////////////////////////////////////////////////
void resize(int w, int h)
{
	GLOBAL_width = w;
	GLOBAL_height = h;
 
	glViewport(0, 0, w, h);
	glOrtho(0, w, 0, h, -1, 1);
 
	set_texture();
}

////////////////////////////////////////////////////////////////////////
void mouseclick(int button, int state, int x, int y)
{
	if (state != GLUT_UP) return;
	
	GLOBAL_cx += (x - GLOBAL_width / 2) * GLOBAL_scale;
	GLOBAL_cy -= (y - GLOBAL_height/ 2) * GLOBAL_scale;
 
	switch(button) {
	case GLUT_LEFT_BUTTON: /* zoom in */
		if (GLOBAL_scale > fabs(x) * 1e-16 && GLOBAL_scale > fabs(y) * 1e-16)
			GLOBAL_scale /= 2;
		break;
	case GLUT_RIGHT_BUTTON: /* zoom out */
		GLOBAL_scale *= 2;
		break;
	/* any other button recenters */
	}

	
	set_texture();

	//print_menu(); // uncomment for convenience; comment for benchmarking
}

////////////////////////////////////////////////////////////////////////
void keypress(unsigned char key, int x, int y)
{
	static int zoomin_x=0, zoomin_y=1; // RUF
	// Ruf: where to start fetching mouse coordinates from GLOBAL_zoomin
	//      (first x coordinate is at GLOBAL_zoomin[0])
	//      (first y coordinate is at GLOBAL_zoomin[1])
	//      (next coordinates are at distance 2: see +=2 at 'z' and 'Z" bellow)
	
	switch(key) {
	case 'q': 
             glFinish();
             glutDestroyWindow(GLOBAL_gwin);
             free(GLOBAL_tex);
             return;
	
	case 27: // Esc
	         GLOBAL_scale = 1./256;
	         GLOBAL_cx = -.6;
	         GLOBAL_cy = 0; 
	         break;
 
	case 'r':
             GLOBAL_color_rotate = (GLOBAL_color_rotate + 1) % 6;
              break;

	case '>':
    case '.':
             GLOBAL_max_iter += 128;
             if (GLOBAL_max_iter > 1 << 15) GLOBAL_max_iter = 1 << 15;
             printf("max iter: %d\n", GLOBAL_max_iter);
             break;
 
	case '<':
    case ',':
             GLOBAL_max_iter -= 128;
             if (GLOBAL_max_iter < 128) GLOBAL_max_iter = 128;
             printf("max iter: %d\n", GLOBAL_max_iter);
             break;

	case 'c':
             GLOBAL_saturation = 1 - GLOBAL_saturation;
             break;
 
	case 's':screen_dump(); return;

	case 'I':GLOBAL_max_iter = 4096; break;
	
	case 'i':GLOBAL_max_iter = 256; break;
	
	case ' ':GLOBAL_invert = !GLOBAL_invert; break;
	
	case 'z':// simulate one mouse click in order to dive one time in zoomin
             GLOBAL_refresh=1;
             mouseclick(GLUT_LEFT_BUTTON, GLUT_UP, GLOBAL_zoomin[zoomin_x], GLOBAL_zoomin[zoomin_y]);
             zoomin_x+=2; zoomin_y+=2;
             break;

	case 'Z':// simulate many mouse clicks in order to dive fully in zoomin

			// COMEÇAR A CONTAR O TEMPO
			gettimeofday(&start,NULL);

             GLOBAL_refresh=1; // use 0 to avoid refreshing all but the last one
             for (zoomin_x=0, zoomin_y=1; zoomin_x < GLOBAL_zoomin_num_pairs; zoomin_x+=2, zoomin_y +=2) {
            	if (zoomin_x == GLOBAL_zoomin_num_pairs-2) GLOBAL_refresh=1;
                
				mouseclick(GLUT_LEFT_BUTTON, GLUT_UP, GLOBAL_zoomin[zoomin_x], GLOBAL_zoomin[zoomin_y]);
					
             }


             // simulate case 's'
             keypress('s', -1, -1);

             // simulate case 'q'
             keypress('q', -1, -1);
			 
			 MPI_Finalize();	

             return;
	}
	set_texture();
	print_menu();
}

////////////////////////////////////////////////////////////////////////
void render()
{
	double	x = (double)GLOBAL_width /GLOBAL_tex_w,
		    y = (double)GLOBAL_height/GLOBAL_tex_h;
 
	glClear(GL_COLOR_BUFFER_BIT);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
 
	glBindTexture(GL_TEXTURE_2D, GLOBAL_texture);
 
	glBegin(GL_QUADS);
 
	glTexCoord2f(0, 0); glVertex2i(0, 0);
	glTexCoord2f(x, 0); glVertex2i(GLOBAL_width, 0);
	glTexCoord2f(x, y); glVertex2i(GLOBAL_width, GLOBAL_height);
	glTexCoord2f(0, y); glVertex2i(0, GLOBAL_height);
 
	glEnd();
 
	glFlush();
	glFinish();
}

////////////////////////////////////////////////////////////////////////

void init_gfx(int *c, char **v, int rank)
{
	if(rank == 0){
		glutInit(c, v);
		glutInitDisplayMode(GLUT_RGB);
		glutInitWindowSize(GLOBAL_window_width, GLOBAL_window_height);
	
		GLOBAL_gwin = glutCreateWindow("Mandelbrot");
		glutDisplayFunc(render);
	
		glutKeyboardFunc(keypress);
		glutMouseFunc(mouseclick);
		glutReshapeFunc(resize);
		glGenTextures(1, &GLOBAL_texture);
	}
	set_texture();
}


int main(int c, char **v)
{

	MPI_Init(&c, &v);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 
	if(rank == 0) // Task 0 executa tudo, incluino o calc mandel
	{
		init_gfx(&c, v, rank);
		print_menu();
		glutMainLoop();	
	} else { // Outras tasks executam só o calc mandel
		for(int c=0; c<20; c++)
			calc_mandel();
		
	}

	MPI_Finalize();	
	
	return 0;

}
