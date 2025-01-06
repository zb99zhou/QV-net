# QV-net: Decentralized Self-Tallying Quadratic Voting with Maximal Ballot Secrecy

This repository provides an implementation of QV-net, a decentralized self-tallying quadratic voting with maximal ballot secrecy.

Decentralized e-voting enables secure and transparent elections without relying on trusted authorities, with blockchain emerging as a popular platform. Quadratic voting (QV) is a social choice mechanism where the cost of votes increases quadratically with the number of votes. To cast $n$ votes, a voter needs to spend $n^2$ tokens, which effectively prevents wealthy token holders from dominating the voting outcome while still allowing them to express strong preferences. 
In Decentralized Autonomous Organizations (DAOs), QV has been adopted by organizations such as MetFi DAO, Karmaverse, Pistachio DAO, Gitcoin, and MoonDAO.

QV-net is the first decentralized quadratic voting scheme, in which voters do not need to trust any external party other than themselves for ballot secrecy. QV-net is self-tallying with maximal ballot secrecy. Self-tallying allows anyone to compute election results once all ballots are cast. Maximal ballot secrecy ensures that what each voter learns from QV-net is nothing more than the tally and their own ballot.


## Tests

Run tests using  `cargo test --release`.