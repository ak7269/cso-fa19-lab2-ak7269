/* NYU Computer Systems Organization Lab 2
 * Rabin-Karp Substring Matching
*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <getopt.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <ctype.h>

#include "rkgrep.h"
#include "bloom.h"

#define PRIME 961748941

// calculate modulo addition, i.e. (a+b) % PRIME
long long
madd(long long a, long long b)
{
	return (a + b) % PRIME;
}

// calculate modulo substraction, i.e. (a-b) % PRIME
long long
msub(long long a, long long b)
{
	return (a>b)?(a-b):(a+PRIME-b);
}

// calculate modulo multiplication, i.e. (a*b) % PRIME
long long
mmul(long long a, long long b)
{
	return (a*b) % PRIME;
}

/* naive_substring_match returns number of positions in the document where
 * the pattern has been found.  In addition, it stores the first position 
 * where the pattern is found in the variable pointed to by first_match_ind
 *
 * Both its arguments "pattern" and "doc" are null-terminated C strings.
 */

int
naive_substring_match(const char *pattern, const char *doc, int *first_match_ind)
{
	int c=0;
	int f=0;
	int l=strlen(pattern);
	int len=strlen(doc);
	int s=len-l;
	int j;
	int i;
	for(j=0;j<=s;j++)
	{
		char temp[l];

		for(i=0;i<l;i++)
		{
	 		temp[i]=doc[i+j];//extracting each of the substrings of length l from the doc and storing it in a char array 
	 	}
	 	temp[l]='\0';//as every string ends with a null character
	 	int a=strcmp(temp,pattern);
	 	if(a==0)
	 	{
	 		if(f==0)
	 		{
				f=1;//flag variable 
				*first_match_ind=j;
			}
			c++;
		}
	}
	return c;
}

/* initialize the Rabin-karp hash computation by calculating 
 * and returning the RK hash over a charbuf of m characters, 
 * i.e. The return value should be 
 * 256^(m-1)*charbuf[0] + 256^(m-2)*charbuf[1] + ... + charbuf[m-1],
 * where 256^(m-1) means 256 raised to the power m-1.
 * rkhash_init stores 256^m in the variable pointed to by h.
 *
 * Note that all operations *, +, - are modulo arithematic, so you 
 * should use the provided functions mmul, madd, msub.
 * (We use "long long" to represent an RK hash)
 */

long long 
rkhash_init(const char *charbuf, int m, long long *h)
{

	int k=0;
	long long hash_value=0;
	int a;
	*h=1;
	long long current;
	for(k=0;k<=m-1;k++)
	{
		for(a=1;a<m-k;a++)
		{
			*h=mmul(*h,256);
		}
		current=mmul(*h,charbuf[k]);//finding (each character*256^(m-k))
		hash_value=madd(hash_value,current);//adding the values to calculate the hash value of the string 
		*h=1;
	}
	return hash_value;
}


/* Given the rabin-karp hash value (curr_hash) over substring Y[i],Y[i+1],...,Y[i+m-1]
 * calculate the hash value over Y[i+1],Y[i+2],...,Y[i+m] = curr_hash * 256 - leftmost * h + rightmost
 * where h is 256 raised to the power m (and given as an argument).  
 * Note that *,-,+ refers to modular arithematic so you should use mmul, msub, madd.
 */
long long 
rkhash_next(long long curr_hash, long long h, char leftmost, char rightmost)
{
	curr_hash=mmul(256,curr_hash);
	curr_hash=msub(curr_hash,mmul(leftmost,h));
	curr_hash=madd(curr_hash,rightmost);

	return curr_hash; 
}

/* rk_substring_match returns the number of positions in the document "doc" where
 * the "pattern" has been found, using the Rabin-karp substring matching algorithm.
 * Both pattern and doc are null-terminated C strings. The function also stores
 * the first position where pattern is found in the int variable pointed to by first_match_ind
 *
 * Note: You should implement the Rabin-Karp algorithm by completing the 
 * rkhash_init and rkhash_next functions and then use them here.
*/
int
rk_substring_match(const char *pattern, const char *doc, int *first_match_ind)
{

	int l = strlen(doc);
	int len = strlen(pattern);
	int p=l-len;
	int i; 
	int c=0,f=0;
	long long h=1;
	long long chash = rkhash_init(doc,len,&h);//storing the hash value of the first substring of the doc 
	long long phash = rkhash_init(pattern,len,&h);

	for(i=1;i<len+1;i++)
	{
		h=mmul(h,256);//calculating 256^len
	}

	for(i=0;i<p;i++)
	{
		if(chash==phash)//if equal the counter increases as the hash value of the pattern is the same as the current substring of the doc's hash value 
		{
			if(f==0)
			{
				f=1;
				*first_match_ind=i;
			}
			c++;
		} 
		chash=rkhash_next(chash,h,doc[i],doc[i+len]);
	}
	return c;
}


/* rk_create_doc_bloom returns a pointer to a newly created bloom_filter. 
 * The new bloom filter is populated with all n-m+1 rabin-karp hashes for 
 * all the substrings of length m in "doc".
 * Hint: use the provided bloom_init() and your implementation of bloom_add() here
 */
bloom_filter *
rk_create_doc_bloom(int m, const char *doc, int bloom_size)
{
	int j;
	bloom_filter *b=bloom_init(bloom_size);//initializing the bloom filter before adding the values
	long long h=1;
	long long current=rkhash_init(doc,m,&h);
	for(j=1;j<m+1;j++)
	{
		h=mmul(256,h);//calculating 256^m where m is the pattern length 
	}
	int l=strlen(doc)-m+1;
	for(int j=0;j<l;j++)
	{
		bloom_add(b, current);
		current=rkhash_next(current,h,doc[j],doc[j+m]);
	}
	return b;
}

/* rk_substring_match_using_bloom returns the total number of positions where "pattern" 
 * is found in "doc".  It performs the matching by first checking against the 
 * pre-populated bloom filter "bf" (which has been created by rk_create_doc_bloom on "doc")
 * If the pattern is not found in "bf", then the function immediately returns 0.
 * Otherwise, the function invokes rk_substring_match() function to find "pattern" in "doc".
*/
int
rk_substring_match_using_bloom(const char *pattern, const char *doc, bloom_filter *bf, int *first_match_ind)
{
	int c=0;
	long long h=1;
	int len=strlen(pattern);
	long long PatHash=rkhash_init(pattern,len,&h);//receiving the pattern's hash value 

	if(bloom_query(bf,PatHash)==true)//checking if the value of this pattern is found in the bloom filter 
	{
		c=rk_substring_match(pattern,doc,first_match_ind);//the number of times the pattern is found in the doc
		return c;
	}
	return 0;	
}