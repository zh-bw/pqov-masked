#ifndef MASK_OV_H
#define MASK_OV_H

#include "src/ov_keypair.h"
#include "gadgets.h"
#include <stdint.h>
#include <stdio.h> 

///////////////////////////// Masked Sign ////////////////////////////////
///
/// @brief Masked signing function for classic secret key.
///
/// @param[out] signature - the signature.
/// @param[in]  masked_sk        - the secret key.
/// @param[in]  message   - the message to be signed.
/// @param[in]  mlen      - the length of the message.
/// @return 0 for success. -1 otherwise.
///
int masked_ov_sign( uint8_t *signature, const masked_sk_t *masked_sk, const uint8_t *message, size_t mlen );

////////////////////////////////////


///
/// @brief Generate key pair for pkc VARIANT
///
/// @param[out] pk        - the compressed public key.
/// @param[out] masked_sk        - the masked secret key.
/// @param[in]  masked_sk_seed   - masked seed for generating secret key.
/// @return 0 for success. -1 otherwise.
///
int masked_generate_keypair_pkc( cpk_t *pk, masked_sk_t *masked_sk, const unsigned char *masked_sk_seed );



#endif