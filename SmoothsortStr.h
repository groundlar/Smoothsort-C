/************************************************************
 *
 * File: SmoothsortStr.h
 * Author: Skylar Cook (sky.spammed@gmail.com)
 * API/external functions for Dijkstra's Smoothsort algorithm, used for
 * string arrays.
 *
 * Thanks goes to Keith Schwarz for his writeup at 
 * http://www.keithschwarz.com/smoothsort/
 *
 * */

#ifndef SMOOTHSORT_STR_H
#define SMOOTHSORT_STR_H


/**
 * Function: void smoothsortFile(*FILE inputfile, *FILE outputfile)
 * ---------------------------------------------------------
 * Given file pointers to the desired input and output files, 
 * sorts the strings found in the input file alphabetically and
 * prints the result into the specified output file. Input format
 * should be one string per line.
 */
void smoothsortFile(FILE *inputfile, FILE *outputfile);


/**
 * Function: void smoothsortArray(char *(strArray)[])
 * ---------------------------------------------------------
 * Given an array of strings and its length as a size_t type, sorts the 
 * array in-place alphabetically.
 */
void smoothsortArray(char *strArray[], size_t arrayLen);

#endif
