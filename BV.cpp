#include "RLWE_utils.h"

#define n 2048
#define l 32
#define Q 61
#define T 2030
#define sigma 8

using namespace std;

pari_sp ltop, lbot;

void print(GEN x){
    cout << GENtostr(x) << endl;
    return;
}

class cryptosystem{
public:
    GEN sk, pk;
    GEN q, t, F;
    
    cryptosystem(){
        GEN v = cgetg(n + 2, t_VEC);
        gel(v, 1) = stoi(1);
        for (int i = 2; i <= n; i++)
            gel(v, i) = gzero;
        gel(v, n + 1) = stoi(1);
        F = gtopolyrev(v, -1);
        q = gshift(stoi(1), Q);
        t = stoi(T);
        while(true){
            q = gnextprime(gadd(q, stoi(1)));
            v = lift(gmodulo(q, gmul(stoi(2), stoi(n))));
            if(gcmp(v, stoi(1)) == 0)
                break;
        }
        t = gnextprime(gadd(t, stoi(1)));
        
        
        generate_secret_key();
        generate_public_key();
        
        return;
    }
    
    GEN sample_error_polynomial(){
        ltop = avma;
        GEN tmp = Sample(n, sigma);
        tmp = gmodulo(tmp, q);
        tmp = gtopolyrev(tmp, -1);
        tmp = gmodulo(tmp, F);
        lbot = avma;
        return tmp;
    }
    
    void generate_secret_key(){
        sk = sample_error_polynomial();
        //gerepile(ltop, lbot, NULL);
        return;
    }
    
    void generate_public_key(){
        int i;
        GEN r, tmp, e;
        pk = cgetg(3, t_VEC);
        
        tmp = cgetg(n + 1, t_VEC);
        for (i = 0; i < n; i++) {
            gel(tmp, i + 1) = generate_random(Q);
        }
        tmp = gmodulo(tmp, q);
        r = gtopolyrev(tmp, -1);
        r = gmodulo(r, F);
        gel(pk, 2) = r;
        
        e = sample_error_polynomial();
        //gerepile(ltop, lbot, NULL);
        
        r = gmul(r, sk);
        e = gmul(t, e);
        gel(pk, 1) = gadd(r, e);
        return;
    }
    
    GEN encrypt(GEN m) {
        GEN ct = cgetg(3, t_VEC);
        GEN u, f, g;
        GEN tmp, M;
        
        u = sample_error_polynomial();
        //gerepile(ltop, lbot, NULL);
        f = sample_error_polynomial();
        //gerepile(ltop, lbot, NULL);
        g = sample_error_polynomial();
        //gerepile(ltop, lbot, NULL);
        
        M = gtovecrev(lift(m));
        M = lift(M);
        M = gmodulo(M, q);
        M = gtopolyrev(M, -1);
        M = gmodulo(M, F);
        
        tmp = gmul(compo(pk, 1), u);
        tmp = gadd(gmul(t, g), tmp);
        gel(ct, 1) = gadd(tmp, M);
        tmp = gmul(compo(pk, 2), u);
        tmp = gadd(gmul(t, f), tmp);
        gel(ct, 2) = gneg(tmp);
        return ct;
    }
    
    GEN decrypt(GEN ct) {
        GEN m;
        GEN tmp, tmp1;
        int ct_len = (int) glength(ct);
        tmp = compo(ct, 1);
        for (int i = 1; i < ct_len; i++) {
            tmp1 = gmul(compo(ct, i + 1), powgi(sk, stoi(i)));
            tmp = gadd(tmp, tmp1);
        }
        
        //print(tmp);
        tmp = lift(tmp);
        tmp1 = lift(gtovecrev(tmp));
        for (int i = 0; i < n; i++) {
            tmp = gshift(compo(tmp1, i + 1), 1);
            if (gcmp(tmp, q) == 1) {
                gel(tmp1, i + 1) = gsub(compo(tmp1, i + 1), q);
            }
        }
        m = lift(gmodulo(tmp1, t));
        return m;
    }
    
    GEN addition(GEN ct_1, GEN ct_2){
        GEN ct;
        int i, k = (int) glength(ct_1), p = (int) glength(ct_2);
        
        if (k >= p) {
            ct = cgetg(k + 1, t_VEC);
            for (i = 0; i < p; i++)
                gel(ct, i + 1) = gadd(compo(ct_1, i + 1), compo(ct_2, i + 1));
            for (i = p; i < k; i++)
                gel(ct, i + 1) = compo(ct_1, i + 1);
        } else {
            ct = cgetg(p + 1, t_VEC);
            for (i = 0; i < k; i++)
                gel(ct, i + 1) = gadd(compo(ct_1, i + 1), compo(ct_2, i + 1));
            for (i = k; i < p; i++)
                gel(ct, i + 1) = compo(ct_2, i + 1);
        }
        return ct;
    }
    
    GEN multiplication(GEN ct_1, GEN ct_2){
        GEN ct;
        int k = (int) glength(ct_1), p = (int) glength(ct_2);
        GEN tmp;
        ct = cgetg(k + p, t_VEC);
        bool flag[k + p - 1];
        for(int i = 0; i < k + p - 1; i++)
            flag[i] = false;
        
        for(int i = 0; i < k; i++)
            for(int j = 0; j < p; j++){
                //cout << i << " " << j << " " << i + j << endl;
                tmp = gmul(compo(ct_1, i + 1), compo(ct_2, j + 1));
                if(flag[i + j])
                    gel(ct, i + j + 1) = gadd(compo(ct, i + j + 1), tmp);
                else{
                    gel(ct, i + j + 1) = tmp;
                    flag[i + j] = true;
                }
            }
        return ct;
    }
    
};

class ciphertext{
public:
    GEN val;
    int degree;
    cryptosystem* pkc;
    
    ciphertext(){};
    
    ciphertext(GEN m, cryptosystem* PKC){
        pkc = PKC;
        val = PKC->encrypt(m);
        degree = 2;
    }
    
    void initialize(GEN m, cryptosystem* PKC){
        pkc = PKC;
        val = PKC->encrypt(m);
        degree = 2;
    }
    
    ciphertext operator+(ciphertext &ct){
        ciphertext result;
        result.val = pkc->addition(this->val, ct.val);
        degree = max(this->degree, ct.degree);
        result.pkc = this->pkc;
        return result;
    }
    
    ciphertext operator*(ciphertext &ct){
        ciphertext result;
        result.val = pkc->multiplication(this->val, ct.val);
        degree = this->degree + ct.degree - 1;
        result.pkc = this->pkc;
        return result;
    }

    GEN decrypt(){
        GEN m = pkc->decrypt(val);
        return m;
    }
};

int main(void) {
    pari_init(2000000000, 2);
    srand(time(NULL));
    cryptosystem pkc;
    GEN message_a, message_b;
    message_a = stoi(10);
    message_b = cgetg(n + 1, t_VEC);
    for(int i = 0; i < n; i++)
        gel(message_b, i + 1) = stoi(i + 1);
    ciphertext a(message_a, &pkc);
    //print(a.decrypt());
    ciphertext b(message_b, &pkc);
    //print(b.decrypt());
    ciphertext result = a * b;
    print(result.decrypt());
    pari_close();
    return 0;
}
