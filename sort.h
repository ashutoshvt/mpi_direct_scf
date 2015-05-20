using namespace psi;
class sort_pred{
public:
boost::shared_ptr<BasisSet> aoBasis; 
sort_pred(boost::shared_ptr<BasisSet> aoBasis)
{
this->aoBasis = aoBasis;
}
bool operator() (const std::pair<int,int> &left, const std::pair<int,int> &right)
         {
        return aoBasis->shell(left.first).nfunction() * aoBasis->shell(left.second).nfunction() >  aoBasis->shell(right.first).nfunction() * aoBasis->shell(right.second).nfunction();
        }
};
